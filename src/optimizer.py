import sys
from typing import List, Dict

import numpy as np
from matplotlib import pyplot as plt
from scipy.linalg import block_diag
from math import inf
import mosek

from read_data import CTMsParameters
from model.control import TrafficControlInput


class MyTrafficOptimizer:
    """
    Class to define optimization problem
    """
    parameters: CTMsParameters      # CTM-s Parameters Object
    phi_0: np.ndarray               # Real Flow Data
    x_0: np.ndarray                 # Initial Condition  # TODO: ???
    k_int: int                      # Length of Time Interval # TODO: Units?
    dt: float                       # Time Increment TODO: Units

    n_c: int                        # Number of Cells
    n_on: int                       # Number of On-ramps
    n_off: int                      # Number of Off-ramps
    n_s: int                        # Number of Stations
    n_r: int                        # Number of Service Station Out-Flows
    m: int                          # Number of Cells + Number of Service Station Out-Flows

    num_constraints: int            # [Mosek]: Number of Constraints in LP
    num_variables: int              # [Mosek]: Number of Optimization Variables in LP

    # CTM-s Dynamics Input to State Matrices # TODO: Come back and fix the maps
    st_to_r: np.ndarray             # Map: ind
    r_to_st: List                   # TODO: ...
    j_to_r: Dict                    # Map: j -> r; cell j -> index of input

    beta: np.ndarray                # TODO: ...
    beta_s: np.ndarray              # TODO: ...
    rho_matrix: np.ndarray          # Input to State Matrix (Flows to Cell Density)
    l_matrix: np.ndarray            # Input to State Matrix (Flows to Number of SS Users)
    e_matrix: np.ndarray            # Input to State Matrix (Flows to SS Queue Length)

    # Linear Program Constraints
    a_ub: np.ndarray                # [Mosek]: Constraint Matrix
    asub: List                      # [Mosek]: Non-zero indices of columns in constraint matrix
    aval: List                      # [Mosek]: Non-zero values of columns in constraint matrix
    b_ub: List                      # [Mosek]: Upper Constraint
    bkc: List

    # Cost Function
    c: np.ndarray                   # [Mosek]: TODO: ...

    # Optimization Parameters
    eta: float                      # Tradeoff between TTS and TTD

    # Control
    control: TrafficControlInput    # TODO: ...
    controlled: bool

    # Solution
    u_opt: np.ndarray               # TODO: ...
    u_opt_k: np.ndarray             # TODO: ...

    # Output Variables
    y_rho: np.ndarray               # TODO: ...
    y_l: np.ndarray                 # TODO: ...
    y_e: np.ndarray                 # TODO: ...
    u_phi: np.ndarray               # TODO: ...
    u_rs: np.ndarray                # TODO: ...

    def __init__(self, parameters, data, k_int, eta):
        self.parameters = parameters
        self.phi_0 = data.tolist()
        self.x_0 = 0  # TODO: Change this
        self.k_int = k_int
        self.dt = self.parameters.highway.dt[0]

        self.n_c = len(self.parameters.highway)
        self.n_on = len(self.parameters.onramps)
        self.n_off = len(self.parameters.offramps)
        self.n_s = len(self.parameters.stations)
        self.n_r = len(self.parameters.stations.j_r)
        self.m = self.n_c + self.n_r
        self.x_0 = np.zeros(self.n_c + 2*self.n_r)  # TODO: Change this...

        self.num_constraints = 0
        self.num_variables = k_int * self.m  # TODO: Should also include on-ramps, off-ramps...

        self.st_to_r = np.zeros_like(self.parameters.stations.j)
        self.r_to_st = [[] for _ in self.parameters.stations.j_r]
        self.j_to_r = {}
        self.build_station_outflow_mapping()
        self.beta = np.zeros((self.k_int, self.n_c))
        self.beta_s = self.parameters.stations.beta_s

        self.rho_matrix = np.zeros((k_int * self.n_c, self.num_variables))
        self.l_matrix = np.zeros((k_int * self.n_s, self.num_variables))
        self.e_matrix = np.zeros((k_int * self.n_s, self.num_variables))

        self.a_ub = []
        self.asub = []
        self.aval = []
        self.b_ub = []
        self.bkc = []

        self.c = np.zeros(self.num_variables)

        self.eta = eta

        self.control = TrafficControlInput()
        self.controlled = False

        self.u_opt = np.zeros(self.num_variables)
        self.u_opt_k = np.zeros((self.k_int, self.m))

        self.y_rho = np.zeros((self.k_int, self.n_c))
        self.y_l = np.zeros((self.k_int, self.n_s))
        self.y_e = np.zeros((self.k_int, self.n_s))
        self.u_phi = np.zeros((self.k_int, self.n_c))
        self.u_rs = np.zeros((self.k_int, self.n_r))

        self.compute_total_beta()
        self.compute_state_matrices()
        self.define_constraints()
        self.define_cost_function()

    def build_station_outflow_mapping(self):
        """
        TODO: Add description / comments
        """
        for id, j in enumerate(self.parameters.stations.j):
            for r_, j_r in enumerate(self.parameters.stations.j_r):
                if j == j_r:
                    self.st_to_r[id] = r_
                    self.r_to_st[r_].append(id)
                    self.j_to_r[j] = np.array([r_])
                    break

    def compute_total_beta(self):
        """
        TODO: Add description / comments
        """
        # TODO: Add comment
        if len(self.parameters.stations.beta_s.shape) == 1:
            k_ones = np.ones((self.k_int, 1))
            self.beta_s = np.kron(k_ones, self.parameters.stations.beta_s)

        for k in range(self.k_int):
            for i_, i in enumerate(self.parameters.stations.i):
                self.beta[k, i] += self.parameters.stations.beta_s[i_]

            for i_, i in enumerate(self.parameters.offramps.i):
                self.beta[k, i] += self.parameters.offramps.beta_r[i_]

    def compute_state_matrices(self):
        """
        TODO: Add description / comments
        """

        a_i = np.divide(self.dt, self.parameters.highway.l)
        a_i_k = np.kron(np.ones((self.k_int, 1)), a_i)
        b_i = np.divide(a_i_k, (1-self.beta))

        #############################
        # Cell Density (Rho Matrix) #
        #############################

        rho_1 = np.diag(a_i)
        rho_4 = np.zeros((self.n_c, self.n_r))

        for k in range(1, self.k_int):
            rho_2 = np.diag(b_i[k-1, :-1], 1)
            rho_3 = rho_1 - rho_2
            rho_5 = np.concatenate((rho_3, rho_4), axis=1)

            for i_, j_r in enumerate(self.parameters.stations.j_r):
                c = self.n_c + self.st_to_r[i_]
                rho_5[j_r, c] += a_i[j_r]

            rho_6 = np.ones((self.k_int-k, 1))
            rho = np.kron(rho_6, rho_5)
            self.rho_matrix[k*self.n_c:, (k-1)*self.m:k*self.m] = rho

        ###############################################
        # Service Station Number of People (L Matrix) #
        ###############################################

        b_i_s = (self.dt * self.beta_s) / (1-self.beta[:, self.parameters.stations.i])

        for k in range(1, self.k_int):
            l_1 = np.zeros((self.n_s, self.m))
            for i_, i in enumerate(self.parameters.stations.i):
                l_1[i_, i] += b_i_s[k-1, i_]
                l_1[i_, self.n_c+self.st_to_r[i_]] -= self.dt

            l_2 = np.ones((self.k_int-k, 1))
            l = np.kron(l_2, l_1)
            self.l_matrix[k*self.n_s:, (k-1)*self.m: k*self.m] = l

        ###########################################
        # Service Station Queue Length (E Matrix) #
        ###########################################

        for k in range(1, self.k_int):
            e_1 = np.zeros((self.n_s, self.m))

            for i_ in self.parameters.stations.id:
                e_1[i_, self.n_c + self.st_to_r[i_]] -= self.dt

            e_2 = np.ones((self.k_int - k, 1))
            e = np.kron(e_2, e_1)
            self.e_matrix[k*self.n_s:, (k-1) * self.m: k*self.m] += e

            for i_, dq in enumerate(self.parameters.stations.delta):
                if k >= dq:
                    for k_dq in range(0, k+1-round(dq)):
                        self.e_matrix[(k*self.n_s) + i_, :] += self.l_matrix[(k_dq*self.n_s) + i_, :]

    def add_constraints(self, a: np.ndarray, b: np.ndarray, c=0):
        assert (a.shape[0] == b.shape[0])

        self.a_ub += a.tolist()
        self.b_ub += b.tolist()

        if c:
            self.bkc += [mosek.boundkey.fx for _ in range(len(b))]
        else:
            self.bkc += [mosek.boundkey.up for _ in range(len(b))]

    def define_constraints(self):
        """
        TODO: Add description / comments
        """
        k_ones = np.ones((self.k_int, 1))

        ############################
        # Queue Length Constraints #
        ############################

        # TODO: Add comments, test
        a_ub_0 = self.e_matrix
        b_ub_0 = np.kron(k_ones, self.parameters.stations.e_max.reshape(-1, 1))

        self.add_constraints(a_ub_0, b_ub_0)

        ###############################
        # Mainstream Flow Constraints #
        ###############################

        # 1) Less than previous demand
        # a) Less than max flow supported from free-flow in previous cell
        beta_ms_i = 1 - self.beta[:, :-1]
        beta_ms_i_k = beta_ms_i.reshape(-1, 1)
        v_i = self.parameters.highway.v[:-1]
        v_i_k = np.kron(k_ones, np.reshape(v_i, (-1, 1)))

        a_1_1 = np.zeros((self.n_c-1, self.m))
        a_1_1[:, 1:self.n_c] = np.eye(self.n_c-1)
        a_1_a_lhs = block_diag(*[a_1_1 for _ in range(self.k_int)])
        c = [self.n_c*k + i for k in range(self.k_int) for i in range(0, self.n_c-1)]
        a_1_a_rhs = np.multiply(np.multiply(beta_ms_i_k, v_i_k), self.rho_matrix[c, :])

        a_ub_1_a = a_1_a_lhs - a_1_a_rhs
        b_ub_1_a = np.zeros(((self.n_c-1) * self.k_int, 1))

        # # Update with initial condition (Densities at k = 0)
        # rho_0 = self.x_0[:self.n_c - 1]
        # b_ub_1_a_0 = np.reshape(np.multiply(np.multiply(beta_ms_i[0, :], v_i), rho_0), (-1, 1))
        # b_ub_1_a[:self.n_c-1, :] = b_ub_1_a_0

        self.add_constraints(a_ub_1_a, b_ub_1_a)

        # b) Less than max flow supported by the previous cell
        q_max_i = self.parameters.highway.q_max[:-1]

        a_ub_1_b = a_1_a_lhs
        b_ub_1_b = np.kron(k_ones, q_max_i.reshape(-1, 1))
        self.add_constraints(a_ub_1_b, b_ub_1_b)

        # 3) Less than priority times supply
        # a) Less than max flow supported due to congestion in next cell
        p_ms_i = self.parameters.highway.p_ms[:-1]
        p_ms_i_k = np.kron(k_ones, p_ms_i.reshape(-1, 1))
        w_i = self.parameters.highway.w[1:]
        w_i_k = np.kron(k_ones, w_i.reshape(-1, 1))
        rho_max_i = self.parameters.highway.rho_max[1:]
        rho_max_i_k = np.kron(k_ones, rho_max_i.reshape(-1, 1))

        a_3_1 = np.zeros((self.n_c - 1, self.m))
        a_3_1[:, 0:self.n_c-1] = np.eye(self.n_c - 1)
        a_3_a_lhs = block_diag(*[a_3_1 for _ in range(self.k_int)])
        a_3_2 = np.multiply(p_ms_i_k, w_i_k)
        c = [self.n_c*k + i for k in range(self.k_int) for i in range(1, self.n_c)]
        a_3_a_rhs = np.multiply(a_3_2, self.rho_matrix[c, :])

        a_ub_3_a = a_3_a_lhs + a_3_a_rhs
        b_ub_3_a = np.multiply(a_3_2, rho_max_i_k)

        # # Update with initial condition (Densities at k = 0)
        # rho_0 = self.x_0[1:self.n_c]
        # b_ub_3_a_0 = np.reshape(np.multiply(w_i, (rho_max_i - rho_0)), (-1, 1))
        # b_ub_3_a[:self.n_c-1, :] = b_ub_3_a_0

        self.add_constraints(a_ub_3_a, b_ub_3_a)

        # b) Less than max flow supported by the next cell
        q_max_i = self.parameters.highway.q_max[1:]
        q_max_i_k = np.kron(k_ones, q_max_i.reshape(-1, 1))

        a_ub_3_b = a_3_a_lhs
        b_ub_3_b = np.multiply(p_ms_i_k, q_max_i_k)
        self.add_constraints(a_ub_3_b, b_ub_3_b)

        ####################################
        # Service Station Flow Constraints #
        ####################################

        # 4) Less than previous demand
        # a) Less than max flow possible from service station
        a_4_1 = np.zeros((self.n_r, self.m))
        a_4_1[:, self.n_c:] = np.eye(self.n_r)
        a_4_a_lhs = block_diag(*[a_4_1 for _ in range(self.k_int)])
        a_4_a_rhs = np.zeros_like(a_4_a_lhs)

        for k in range(self.k_int):
            for r_, stations in enumerate(self.r_to_st):
                c_1 = (k * self.n_r) + r_
                for id in stations:
                    c_2 = (k * self.n_s) + id
                    a_4_a_rhs[c_1, :] += self.e_matrix[c_2, :] / self.dt
                    dq = self.parameters.stations.delta[id]
                    if k >= dq:
                        k_dq = k-round(dq)
                        c_3 = (k_dq * self.n_s) + id
                        a_4_a_rhs[c_1, :] += self.l_matrix[c_3, :] / self.dt

        a_ub_4_a = a_4_a_lhs - a_4_a_rhs
        b_ub_4_a = np.kron(k_ones, np.zeros((self.n_r, 1)))
        self.add_constraints(a_ub_4_a, b_ub_4_a)

        # b) Less than max flow supported by the service station
        ids = [st[0] for st in self.r_to_st]
        r_s_max_i = self.parameters.stations.r_s_max[ids]

        a_ub_4_b = a_4_a_lhs
        b_ub_4_b = np.kron(k_ones, r_s_max_i.reshape(-1, 1))
        self.add_constraints(a_ub_4_b, b_ub_4_b)

        # c) Less than max flow prescribed by controller
        if self.controlled:
            a_ub_4_c = a_ub_4_b
            b_ub_4_c = self.control.r_s_c.reshape(-1, 1)
            self.add_constraints(a_ub_4_c, b_ub_4_c)

        # 6) Less than priority times supply
        # a) Less than max flow supported due to congestion in next cell
        # TODO: Add comments here...
        a_6_1 = self.parameters.stations.j_r + 1
        j_r = a_6_1[a_6_1 < len(self.parameters.highway)] - 1

        if j_r.any():
            p_s_i = np.array([self.parameters.stations.j_to_p[j] for j in j_r])
            p_s_i_k = np.kron(k_ones, p_s_i.reshape(-1, 1))
            w_i = self.parameters.highway.w[j_r+1]
            w_i_k = np.kron(k_ones, w_i.reshape(-1, 1))
            rho_max_i = self.parameters.highway.rho_max[j_r+1]
            rho_max_i_k = np.kron(k_ones, rho_max_i.reshape(-1, 1))

            a_6_3 = np.zeros((len(j_r), self.m))

            for i_, j in enumerate(j_r):
                a_6_3[i_, self.n_c + self.j_to_r[j]] += 1

            a_6_a_lhs = block_diag(*[a_6_3 for _ in range(self.k_int)])

            c = np.array([(self.n_c*k) + (j+1) for k in range(self.k_int) for j in j_r])
            a_6_4 = np.multiply(w_i_k, self.rho_matrix[c, :])
            a_6_a_rhs = np.multiply(p_s_i_k, a_6_4)

            a_ub_6_a = a_6_a_lhs + a_6_a_rhs
            b_ub_6_a = np.multiply(p_s_i_k, np.multiply(w_i_k, rho_max_i_k))

            # # TODO: Update with initial condition (Densities at k = 0)
            # rho_0 = self.x_0[j_r]
            # b_ub_6_a_0 = np.reshape(np.multiply(w_i, (rho_max_i - rho_0)), (-1, 1))
            # b_ub_6_a[:self.n_c - 1, :] = b_ub_6_a_0

            self.add_constraints(a_ub_6_a, b_ub_6_a)

            # b) Less than max flow supported by the next cell
            q_max_i = self.parameters.highway.q_max[j_r + 1]
            q_max_i_k = np.multiply(k_ones, q_max_i.reshape(-1, 1))

            a_ub_6_b = a_6_a_lhs
            b_ub_6_b = np.multiply(p_s_i_k, q_max_i_k)
            self.add_constraints(a_ub_6_b, b_ub_6_b)

            # #################################
            #
            # 5) Less than supply minus mainstream demand
            a_5_lhs = a_6_a_lhs
            a_5_rhs_1 = np.multiply(w_i_k, self.rho_matrix[c, :])
            b_ub_5_1 = np.multiply(p_s_i_k, np.multiply(w_i_k, rho_max_i_k))
            b_ub_5_2 = np.multiply(p_s_i_k, q_max_i_k)
            beta_ms_j = 1 - self.beta[:, j_r]
            beta_ms_j_k = beta_ms_j.reshape(-1, 1)
            v_j = self.parameters.highway.v[j_r]
            v_j_k = np.kron(k_ones, np.reshape(v_j, (-1, 1)))
            a_5_1 = np.multiply(beta_ms_j_k, v_j_k)
            a_5_rhs_3 = np.multiply(a_5_1, self.rho_matrix[c-1, :])
            b_ub_5_3 = np.zeros((self.k_int*len(j_r), 1))
            b_ub_5_4 = np.kron(k_ones, self.parameters.highway.q_max[j_r])

            a_ub_5_a = a_5_lhs + a_5_rhs_1 + a_5_rhs_3
            b_ub_5_a = b_ub_5_1 - b_ub_5_3
            self.add_constraints(a_ub_5_a, b_ub_5_a)

            # a_ub_5_b = a_5_lhs + a_5_rhs_1
            # b_ub_5_b = b_ub_5_1 - b_ub_5_4
            # self.add_constraints(a_ub_5_b, b_ub_5_b)

            a_ub_5_c = a_5_lhs + a_5_rhs_3
            b_ub_5_c = b_ub_5_2 - b_ub_5_3
            self.add_constraints(a_ub_5_c, b_ub_5_c)
            #
            # # a_ub_5_d = a_5_lhs + a_5_rhs_3
            # # b_ub_5_d = b_ub_5_2 - b_ub_5_4
            # # self.add_constraints(a_ub_5_d, b_ub_5_d)
            #
            # 2) Less than supply minus service-station demand
            a_2_1 = np.zeros((len(j_r), self.m))

            for i, j in enumerate(j_r):
                a_2_1[i, j] += 1

            a_2_lhs = block_diag(*[a_2_1 for _ in range(self.k_int)])

            a_2_rhs_1 = a_5_rhs_1
            b_ub_2_1 = b_ub_5_1
            b_ub_2_2 = b_ub_5_2

            a_2_rhs_3 = np.zeros((len(j_r)*self.k_int, self.num_variables))

            for j in j_r:
                r = self.j_to_r[j]
                c_1 = [k*len(j_r) + r[0] for k in range(self.k_int)]
                c_2 = [k*self.n_r + r[0] for k in range(self.k_int)]
                a_2_rhs_3[c_1, :] = a_4_a_rhs[c_2, :]

            b_ub_2_3 = np.kron(k_ones, np.zeros((len(j_r), 1)))
            r = [self.j_to_r[j] for j in j_r]
            ids = [self.r_to_st[r_] for r_ in r[0]]
            r_s_max_j = self.parameters.stations.r_s_max[tuple(ids)]
            b_ub_2_4 = np.kron(k_ones, r_s_max_j)

            a_ub_2_a = a_2_lhs + a_2_rhs_1 + a_2_rhs_3
            b_ub_2_a = b_ub_2_1 - b_ub_2_3
            self.add_constraints(a_ub_2_a, b_ub_2_a)

            # a_ub_2_b = a_2_lhs + a_2_rhs_1
            # b_ub_2_b = b_ub_2_1 - b_ub_2_4
            # self.add_constraints(a_ub_2_b, b_ub_2_b)

            a_ub_2_c = a_2_lhs + a_2_rhs_3
            b_ub_2_c = b_ub_2_2 - b_ub_2_3
            self.add_constraints(a_ub_2_c, b_ub_2_c)

            # a_ub_2_d = a_2_lhs + a_2_rhs_3
            # b_ub_2_d = b_ub_2_2 - b_ub_2_4
            # self.add_constraints(a_ub_2_d, b_ub_2_d)

            if self.controlled:
                b_ub_2_5 = self.control.r_s_c[:, self.j_to_r[j_r]].reshape(-1, 1)
                a_ub_2_e = a_2_lhs + a_2_rhs_1
                b_ub_2_e = b_ub_2_1 - b_ub_2_5
                a_ub_2_f = a_2_lhs
                b_ub_2_f = b_ub_2_2 - b_ub_2_5

                self.add_constraints(a_ub_2_e, b_ub_2_e)
                self.add_constraints(a_ub_2_f, b_ub_2_f)

            # #################################

        # # Positive State Constraints
        # a_ub_7_a = -self.rho_matrix
        # b_ub_7_a = np.zeros((self.n_c*self.k_int, 1))
        # a_ub_7_b = -self.l_matrix
        # b_ub_7_b = np.zeros((self.n_s*self.k_int, 1))
        # a_ub_7_c = -self.e_matrix
        # b_ub_7_c = np.zeros((self.n_s*self.k_int, 1))
        #
        # self.a_ub += a_ub_7_a.tolist()
        # self.b_ub += b_ub_7_a.tolist()
        # self.bkc += [mosek.boundkey.up for _ in range(len(b_ub_7_a))]
        #
        # self.a_ub += a_ub_7_b.tolist()
        # self.b_ub += b_ub_7_b.tolist()
        # self.bkc += [mosek.boundkey.up for _ in range(len(b_ub_7_b))]
        #
        # self.a_ub += a_ub_7_c.tolist()
        # self.b_ub += b_ub_7_c.tolist()
        # self.bkc += [mosek.boundkey.up for _ in range(len(b_ub_7_a))]

        # # Enforce that density of the last cell is equal to 0
        # i = [(k*self.n_c)-1 for k in range(self.k_int)]
        # a_eq = self.rho_matrix[i, :]
        # b_eq = np.zeros((self.k_int, 1))
        #
        # self.a_ub += a_eq.tolist()
        # self.b_ub += b_eq.tolist()
        # self.bkc += [mosek.boundkey.fx for _ in range(len(b_eq))]

        self.formulate_constraints()

    # TODO: Test this function!
    def formulate_constraints(self):
        """
        TODO: Add description / comments
        """

        self.num_constraints = len(self.a_ub)
        # self.num_constraints = len(self.a_ub) + len(self.a_eq)
        self.a_ub = np.array(self.a_ub)

        for column in self.a_ub.T:
            ind = np.nonzero(column)[0].tolist()
            self.asub += [ind]
            self.aval += [column[ind].tolist()]

    # TODO: Test this function!!
    def define_cost_function(self):
        """
        TODO: Add description / comments
        """

        # Define coefficients from Total Time Spent (TTS: Minimize)
        k_ones = np.ones((self.k_int, 1))
        l_i_k = np.reshape(np.kron(k_ones, self.parameters.highway.l), (-1, 1))
        tts_1 = np.sum(np.multiply(l_i_k, self.rho_matrix), axis=0)
        tts_2 = np.sum(self.e_matrix, axis=0)
        tts = self.dt * (tts_1 + tts_2)

        # Define coefficients from Total Travel Distance (TTD: Maximize)
        l_1 = np.concatenate((self.parameters.highway.l, self.parameters.stations.l_r))
        l_2 = np.ones(self.k_int)
        ttd = self.dt * np.kron(l_2, l_1)

        # self.c = tts - (self.eta * ttd)
        self.c = -np.ones(self.num_variables)

    def update_control(self, new_control: TrafficControlInput):
        """
        TODO: Add description / comments
        """

        self.control = new_control
        self.controlled = True

        # TODO: self.compute_total_beta()
        # TODO: self.define_constraints()
        # TODO: self.formulate_constraints()
        # TODO: self.define_cost_function()

    def update_variables(self):
        """
        TODO: Add comments
        """

        self.u_opt_k = np.reshape(self.u_opt, (self.k_int, self.m))

        self.y_rho = np.reshape(np.matmul(self.rho_matrix, self.u_opt), (self.k_int, self.n_c))
        self.y_l = np.reshape(np.matmul(self.l_matrix, self.u_opt), (self.k_int, self.n_s))
        self.y_e = np.reshape(np.matmul(self.e_matrix, self.u_opt), (self.k_int, self.n_s))

        self.u_phi = self.u_opt_k[:, :self.n_c]
        self.u_rs = self.u_opt_k[:, self.n_c:]

    def get_output(self) -> np.ndarray:
        """
        TODO: ...
        """

        return self.y_rho

    # Define a stream printer to grab output from MOSEK
    def streamprinter(self, text):
        """
        TODO: Add description / comments
        """

        sys.stdout.write(text)
        sys.stdout.flush()

    def plot_output(self, cell_ids: List, loc: str):
        """
        TODO: Plot density [veh/km] of cell i during the time interval
        """

        fig = plt.figure()
        t = np.arange(0, self.k_int)
        rho_i = self.y_rho[:, cell_ids]
        plt.plot(t, rho_i)
        fig.savefig(loc + 'lp/rho.png', dpi=300)  # TODO: Naming scheme for files...

    def plot_fundamental_diagram(self, cell_id: int, loc: str):
        """
        TODO: ...
        """
        raise NotImplementedError

        fig = plt.figure()
        t = np.arange(0, self.k_int)
        rho_i = self.y_rho[:, cell_id]
        phi_i = self.u_opt_k[:, cell_id]
        fig.savefig(loc + 'lp/fundamental_diagram.png', dpi=300)  # TODO: Naming scheme for files...

    def solve(self, print_sol: bool):
        """
        TODO: Add description / comments
        """

        # Create a task object
        with mosek.Task() as task:
            # Attach a log stream printer to the task
            task.set_Stream(mosek.streamtype.log, self.streamprinter)

            # Append 'numConstraints' empty constraints.
            # The constraints will initially have no bounds.
            task.appendcons(self.num_constraints)

            # Append 'numVariables' variables.
            # The variables will initially be fixed at zero (x=0).
            task.appendvars(self.num_variables)

            for j in range(self.num_variables):
                # Set the linear term c_j in the objective.
                task.putcj(j, self.c[j])

                # Set the bounds on variable j
                # blx[j] <= x_j <= bux[j]
                (k, r) = divmod(j, self.m)

                if r == 0:  # Initial flow condition
                    task.putvarbound(j,
                                     mosek.boundkey.fx,     # TODO: ...
                                     self.phi_0[k],         # TODO: ...
                                     self.phi_0[k])         # TODO: ...
                else:
                    task.putvarbound(j,                     # TODO: ...
                                     mosek.boundkey.lo,     # TODO: ...
                                     0,                     # TODO: ...
                                     +inf)                  # TODO: ...

                # Input column j of A
                task.putacol(j,                             # Variable (column) index.
                             self.asub[j],                  # Row index of non-zeros in column j.
                             self.aval[j])                  # Non-zero Values of column j.

            # Set the bounds on constraints.
            # blc[i] <= constraint_i <= buc[i]
            for i in range(self.num_constraints):
                task.putconbound(i,                         # TODO: ...
                                 self.bkc[i],               # TODO: ...
                                 self.b_ub[i][0],           # Lower Bound is ignored
                                 self.b_ub[i][0])           # TODO: ...

            # Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.minimize)

            # Solve the problem
            task.optimize()

            # Print a summary containing information
            # about the solution for debugging purposes
            task.solutionsummary(mosek.streamtype.msg)

            # Get status information about the solution
            solsta = task.getsolsta(mosek.soltype.bas)

            if solsta == mosek.solsta.optimal:
                self.u_opt = task.getxx(mosek.soltype.bas)
                self.update_variables()

                if print_sol:
                    print("Optimal solution: ")
                    for i in range(self.num_variables):
                        print("x[" + str(i) + "]=" + str(self.u_opt[i]))

            elif (solsta == mosek.solsta.dual_infeas_cer or
                  solsta == mosek.solsta.prim_infeas_cer):
                print("Primal or dual infeasibility certificate found.\n")
            elif solsta == mosek.solsta.unknown:
                print("Unknown solution status")
            else:
                print("Other solution status")
