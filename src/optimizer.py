import sys
from typing import List, Dict

import numpy as np
from scipy.linalg import block_diag
from math import inf
import mosek

from model.parameters import CTMsParameters
from model.control import TrafficControlInput


class MyTrafficOptimizer:
    """
    Class to define optimization problem
    """
    parameters: CTMsParameters  # CTM-s Parameters Object
    x_0: np.ndarray  # Initial Highway State at time of LP initialization
    phi_0: np.ndarray  # Incoming Flow at First Cell
    phi_s_0: np.ndarray   # TODO: ...

    k_0: int  # Time at LP initialization # TODO?
    k_l: int  # Length of Time Interval # TODO: Units?
    dt: float  # Time Increment TODO: Units

    n_c: int  # Number of Cells
    n_on: int  # Number of On-ramps
    n_off: int  # Number of Off-ramps
    n_s: int  # Number of Stations
    n_r: int  # Number of Service Station Out-Flows
    n: int  # Dimension of State: Density of each cell, users and queue length for each station
    m: int  # Number of flows between cells, exiting service stations, terminal

    num_constraints: int  # [Mosek]: Number of Constraints in LP
    num_variables: int  # [Mosek]: Number of Optimization Variables in LP

    st_to_r: np.ndarray  # CTM-s Input to State Matrices # TODO: Come back and fix the maps, comments
    r_to_st: List  # TODO: ...
    j_to_r: Dict  # Map: j -> r; cell j -> index of input

    beta: np.ndarray  # TODO: ...
    beta_s: np.ndarray  # TODO: ...

    rho_0_matrix: np.ndarray  # TODO: Selection Matrix
    l_0_matrix: np.ndarray  # TODO: Selection Matrix
    l_s_0_matrix: np.ndarray # TODO: ...?
    e_0_matrix: np.ndarray  # TODO: Selection Matrix
    e_s_0_matrix: np.ndarray  # TODO: ???
    rho_matrix: np.ndarray  # Input to State Matrix (Flows to Cell Density)
    l_matrix: np.ndarray  # Input to State Matrix (Flows to Number of SS Users)
    e_matrix: np.ndarray  # Input to State Matrix (Flows to SS Queue Length)

    a_ub: np.ndarray  # [Mosek]: Constraint Matrix
    asub: List  # [Mosek]: Non-zero indices of columns in constraint matrix
    aval: List  # [Mosek]: Non-zero values of columns in constraint matrix
    b_lb: List  # [Mosek]: Lower Constraint
    b_ub: List  # [Mosek]: Upper Constraint
    bkc: List  # [Mosek]: Constraint Bound Keys
    c: np.ndarray  # [Mosek]: Cost Function

    # Optimization Parameters
    eta: float  # Tradeoff between TTS and TTD

    # Control
    control: TrafficControlInput  # TODO: ...
    controlled: bool

    # Solution
    u_opt: np.ndarray  # TODO: ...
    u_opt_k: np.ndarray  # TODO: ...

    # State, Input Variables
    x_rho: np.ndarray  # TODO: ...
    x_l: np.ndarray  # TODO: ...
    x_e: np.ndarray  # TODO: ...
    u_phi: np.ndarray  # TODO: ...
    u_rs: np.ndarray  # TODO: ...

    # Helper
    k_ones: np.ndarray
    k_1_ones: np.ndarray

    def __init__(self, parameters: CTMsParameters, eta: float):
        self.parameters = parameters
        self.x_0 = np.empty(0)
        self.phi_0 = np.empty(0)
        self.phi_s_0 = np.empty(0)
        self.k_0 = 0
        self.k_l = 0
        self.dt = parameters.highway.dt

        self.n_c = len(self.parameters.highway)
        self.n_on = len(self.parameters.onramps)
        self.n_off = len(self.parameters.offramps)
        self.n_s = len(self.parameters.stations)
        self.n_r = len(self.parameters.stations.j_r)
        self.n = self.n_c + 2 * self.n_s
        self.m = self.n_c + self.n_r + 1  # TODO: Include on-ramps eventually

        self.num_constraints = 0
        self.num_variables = self.k_l * self.m

        self.st_to_r = np.zeros_like(self.parameters.stations.j)
        self.r_to_st = [[] for _ in self.parameters.stations.j_r]
        self.j_to_r = {}
        self.build_station_outflow_mapping()

        self.beta = np.zeros((self.k_l, self.n_c))
        self.beta_s = self.parameters.stations.beta_s

        self.rho_0_matrix = np.empty(0)
        self.l_0_matrix = np.empty(0)
        self.l_s_0_matrix = np.empty(0)
        self.e_0_matrix = np.empty(0)
        self.e_s_0_matrix = np.empty(0)

        self.rho_matrix = np.empty(0)
        self.l_matrix = np.empty(0)
        self.e_matrix = np.empty(0)

        self.a_ub = []
        self.asub = []
        self.aval = []
        self.b_lb = []
        self.b_ub = []
        self.bkc = []

        self.c = np.empty(0)

        self.eta = eta

        # TODO: ...
        self.control = TrafficControlInput()
        self.controlled = False

        self.u_opt = np.empty(0)
        self.u_opt_k = np.empty(0)

        self.x_rho = np.empty(0)
        self.x_l = np.empty(0)
        self.x_e = np.empty(0)
        self.u_phi = np.empty(0)  # Including terminal flow
        self.u_rs = np.empty(0)

        self.k_ones = np.empty(0)
        self.k_1_ones = np.empty(0)

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
        # TODO: Time dependent beta?
        if self.parameters.stations.beta_s.ndim == 1:
            self.beta_s = np.kron(self.k_ones, self.parameters.stations.beta_s)

        for k in range(self.k_l):
            for i_, i in enumerate(self.parameters.stations.i):
                self.beta[k, i] += self.beta_s[k, i_]

            for i_, i in enumerate(self.parameters.offramps.i):
                self.beta[k, i] += self.parameters.offramps.beta_r[i_]

    def compute_state_matrices(self):
        """
        TODO: Add description / comments
        """

        a_i = np.divide(self.dt, self.parameters.highway.l)
        a_i_k = np.kron(self.k_ones, a_i)
        b_i = np.divide(a_i_k, (1 - self.beta))

        #############################
        # Cell Density (Rho Matrix) #
        #############################

        self.rho_0_matrix[:, :self.n_c] = np.kron(self.k_1_ones, np.eye(self.n_c))

        rho_1 = np.diag(a_i)
        rho_4 = np.zeros((self.n_c, self.n_r))

        for k in range(1, self.k_l + 1):
            rho_2 = np.diag(b_i[k - 1, :])
            rho_3 = np.zeros((self.n_c, self.n_c + 1))
            rho_3[:, :self.n_c] += rho_1
            rho_3[:, 1:self.n_c + 1] -= rho_2
            rho_5 = np.concatenate((rho_3, rho_4), axis=1)

            for i_, j_r in enumerate(self.parameters.stations.j_r):
                rho_5[j_r, self.n_c + 1 + self.st_to_r[i_]] += a_i[j_r]

            rho_6 = np.ones((self.k_l - k + 1, 1))
            rho = np.kron(rho_6, rho_5)
            self.rho_matrix[k * self.n_c:, (k - 1) * self.m:k * self.m] = rho

        ###############################################
        # Service Station Number of People (L Matrix) #
        ###############################################

        l_0 = np.zeros((self.n_s, self.n))
        for i in self.parameters.stations.id:
            l_0[i, self.n_c + (2 * i)] = 1

        self.l_0_matrix = np.kron(self.k_1_ones, l_0)

        b_i_s = (self.dt * self.beta_s) / (1 - self.beta[:, self.parameters.stations.i])
        b_i_s_0 = b_i_s[0]  # TODO: Should be the time history of beta, come back and fix this...

        for k in range(1, self.k_l + 1):
            l_1 = np.zeros((self.n_s, self.m))

            for i_, i in enumerate(self.parameters.stations.i):
                l_1[i_, i + 1] += b_i_s[k - 1, i_]

            l_2 = np.ones((self.k_l - k + 1, 1))
            l_3 = np.kron(l_2, l_1)
            self.l_matrix[k * self.n_s:, (k - 1) * self.m:k * self.m] += l_3

        for k in range(1, self.k_l + 1):
            for i_, dq in enumerate(self.parameters.stations.delta):
                i = self.parameters.stations.i[i_]
                r_1 = (k * self.n_s) + i_

                k_dq_1 = self.k_0 - round(dq)
                k_dq_2 = k_dq_1 + k

                if k_dq_2 > 0:
                    for k_dq in range(max(k_dq_1, 0), min(k_dq_2, self.k_0)):
                        self.l_s_0_matrix[r_1, k_dq] -= b_i_s_0  # TODO: Change for non-constant beta

                    if k_dq_2 > self.k_0:
                        r_2 = ((k_dq_2 - self.k_0) * self.n_s) + i_
                        self.l_matrix[r_1, :] -= self.l_matrix[r_2, :]

        ###########################################
        # Service Station Queue Length (E Matrix) #
        ###########################################

        e_0 = np.zeros((self.n_s, self.n))
        for i in self.parameters.stations.id:
            e_0[i, self.n_c + 1 + (2 * i)] += 1

        self.e_0_matrix = np.kron(self.k_1_ones, e_0)

        for k in range(1, self.k_l + 1):
            e_1 = np.zeros((self.n_s, self.m))
            for i in self.parameters.stations.id:
                e_1[i, self.n_c + 1 + self.st_to_r[i]] -= self.dt

            e_2 = np.ones((self.k_l - k + 1, 1))
            e_3 = np.kron(e_2, e_1)
            self.e_matrix[k * self.n_s:, (k - 1) * self.m:k * self.m] += e_3

        for k in range(1, self.k_l + 1):
            for i_, dq in enumerate(self.parameters.stations.delta):
                r_1 = (k * self.n_s) + i_

                k_dq_1 = self.k_0 - round(dq)
                k_dq_2 = k_dq_1 + k

                if k_dq_2 > 0:
                    for k_dq in range(max(k_dq_1, 0), min(k_dq_2, self.k_0)):
                        self.e_s_0_matrix[r_1, k_dq] += b_i_s_0    # TODO: Change for non-constant beta
                    if k_dq_2 > self.k_0:
                        r_2 = ((k_dq_2 - self.k_0) * self.n_s) + i_
                        self.e_matrix[r_1, :] += self.l_matrix[r_2, :]

    def add_constraints(self, a: np.ndarray, b_1: np.ndarray, bkc=mosek.boundkey.up, b_2=np.empty(0)):
        """
        TODO: Add description and comments
        """

        assert (a.shape[0] == b_1.shape[0])

        self.a_ub += a.tolist()

        if bkc == mosek.boundkey.ra:
            self.b_lb += b_2.tolist()
            self.b_ub += b_1.tolist()
        else:
            self.b_lb += b_1.tolist()
            self.b_ub += b_1.tolist()

        self.bkc += [bkc for _ in range(len(b_1))]

    def formulate_constraints(self):
        """
        TODO: Add description / comments
        """

        self.num_constraints = len(self.a_ub)
        self.a_ub = np.array(self.a_ub)

        for column in self.a_ub.T:
            ind = np.nonzero(column)[0].tolist()
            self.asub += [ind]
            self.aval += [column[ind].tolist()]

    def define_state_constraints(self):
        """
        TODO: ...
        """

        ################################
        # Max Cell Density Constraints #
        ################################

        rho_max = self.parameters.highway.rho_max.reshape((-1, 1))
        b_rho_0 = np.zeros(((self.k_l + 1) * self.n_c, 1)) - np.matmul(self.rho_0_matrix, self.x_0)
        b_rho_max = np.kron(self.k_1_ones, rho_max) - np.matmul(self.rho_0_matrix, self.x_0)
        self.add_constraints(self.rho_matrix, b_rho_max, mosek.boundkey.ra, b_rho_0)

        ############################
        # Queue Length Constraint #
        ############################

        e_max = self.parameters.stations.e_max.reshape(-1, 1)
        b_e_0 = np.zeros(((self.k_l + 1) * self.n_s, 1)) - np.matmul(self.e_0_matrix, self.x_0) - np.matmul(self.e_s_0_matrix, self.phi_s_0)
        b_e_max = np.kron(self.k_1_ones, e_max) - np.matmul(self.e_0_matrix, self.x_0) - np.matmul(self.e_s_0_matrix, self.phi_s_0)
        self.add_constraints(self.e_matrix, b_e_max, mosek.boundkey.ra, b_e_0)

        ############################################
        # Max Service Station Capacity Constraints #
        ############################################
        b_l_0 = np.zeros(((self.k_l + 1) * self.n_s, 1)) - np.matmul(self.l_0_matrix, self.x_0) - np.matmul(self.l_s_0_matrix, self.phi_s_0)
        self.add_constraints(self.l_matrix, b_l_0, mosek.boundkey.lo)
        # TODO!

    def define_flow_constraints(self):
        """
        TODO: Add description / comments
        """

        ###############################
        # Mainstream Flow Constraints #
        ###############################

        # 1) Less than previous demand
        # a) Less than max flow supported from free-flow in previous cell
        # Note: But no constraint on phi_0, phi_final
        beta_ms_i = 1 - self.beta[:, :-1]
        beta_ms_i_k = beta_ms_i.reshape(-1, 1)
        v_i = self.parameters.highway.v[:-1]
        v_i_k = np.kron(self.k_ones, v_i.reshape(-1, 1))

        a_1_1 = np.zeros((self.n_c - 1, self.m))
        a_1_1[:, 1:self.n_c] = np.eye(self.n_c - 1)
        a_1_a_lhs = block_diag(*[a_1_1 for _ in range(self.k_l)])
        a_1_2 = np.multiply(beta_ms_i_k, v_i_k)
        r_1 = [(self.n_c * k) + i for k in range(self.k_l) for i in range(self.n_c - 1)]
        a_1_a_rhs = np.multiply(a_1_2, self.rho_matrix[r_1, :])

        a_ub_1_a = a_1_a_lhs - a_1_a_rhs
        b_ub_1_a = np.multiply(a_1_2, np.matmul(self.rho_0_matrix[r_1, :], self.x_0))  # Initial condition
        self.add_constraints(a_ub_1_a, b_ub_1_a)

        # b) Less than max flow supported by the previous cell
        # Note: But no constraint on phi_0, phi_final
        q_max_i = self.parameters.highway.q_max[:-1]

        a_ub_1_b = a_1_a_lhs
        b_ub_1_b = np.kron(self.k_ones, q_max_i.reshape(-1, 1))
        self.add_constraints(a_ub_1_b, b_ub_1_b)

        # 3) Less than priority times supply
        # a) Less than max flow supported due to congestion in next cell
        # Note: But no constraint on phi_final
        p_ms_i = self.parameters.highway.p_ms
        p_ms_i_k = np.kron(self.k_ones, p_ms_i.reshape(-1, 1))
        w_i = self.parameters.highway.w
        w_i_k = np.kron(self.k_ones, w_i.reshape(-1, 1))
        rho_max_i = self.parameters.highway.rho_max
        rho_max_i_k = np.kron(self.k_ones, rho_max_i.reshape(-1, 1))

        a_3_1 = np.zeros((self.n_c, self.m))
        a_3_1[:, :self.n_c] = np.eye(self.n_c)
        a_3_a_lhs = block_diag(*[a_3_1 for _ in range(self.k_l)])
        a_3_2 = np.multiply(p_ms_i_k, w_i_k)
        r_3 = [(self.n_c * k) + i for k in range(self.k_l) for i in range(self.n_c)]
        a_3_a_rhs = np.multiply(a_3_2, self.rho_matrix[r_3, :])

        a_ub_3_a = a_3_a_lhs + a_3_a_rhs
        b_ub_3_a = np.multiply(a_3_2, rho_max_i_k - np.matmul(self.rho_0_matrix[r_3, :], self.x_0))  # Initial condition
        self.add_constraints(a_ub_3_a, b_ub_3_a)

        # b) Less than max flow supported by the next cell
        # Note: But no constraint on phi_final
        q_max_i = self.parameters.highway.q_max
        q_max_i_k = np.kron(self.k_ones, q_max_i.reshape(-1, 1))

        a_ub_3_b = a_3_a_lhs
        b_ub_3_b = np.multiply(p_ms_i_k, q_max_i_k)
        self.add_constraints(a_ub_3_b, b_ub_3_b)

        ####################################
        # Service Station Flow Constraints #
        ####################################

        # 4) Less than previous demand
        # a) Less than max flow possible from service station
        a_4_1 = np.zeros((self.n_r, self.m))
        a_4_1[:, self.n_c + 1:] = np.eye(self.n_r)
        a_4_a_lhs = block_diag(*[a_4_1 for _ in range(self.k_l)])
        a_4_a_rhs = np.zeros_like(a_4_a_lhs)

        b_ub_4_a = np.zeros((self.k_l * self.n_r, 1))

        b_i_s = (self.dt * self.beta_s) / (1 - self.beta[:, self.parameters.stations.i])
        b_i_s_0 = b_i_s[0]  # TODO: Change for non-constant beta

        for k in range(self.k_l):
            for r_, stations in enumerate(self.r_to_st):
                r_1 = (k * self.n_r) + r_
                for i_ in stations:
                    r_2 = (k * self.n_s) + i_
                    a_4_a_rhs[r_1, :] += self.e_matrix[r_2, :]
                    b_ub_4_a[r_1] += (np.matmul(self.e_0_matrix[r_2, :], self.x_0) / self.dt)
                    b_ub_4_a[r_1] += np.matmul(self.e_s_0_matrix[r_2, :], self.phi_s_0)

                    dq = self.parameters.stations.delta[i_]
                    k_dq_1 = self.k_0 - round(dq)
                    k_dq_2 = k_dq_1 + k

                    i = self.parameters.stations.i[i_]
                    if k_dq_2 >= 0:
                        if k_dq_2 < self.k_0:
                            b_ub_4_a[r_1] += b_i_s_0 * self.phi_s_0[k_dq_2] / self.dt
                        else:
                            c_1 = ((k_dq_2 - self.k_0) * self.m) + i + 1
                            a_4_a_rhs[r_1, c_1] += (b_i_s[k_dq_2 - self.k_0, i_] / self.dt)

        a_ub_4_a = a_4_a_lhs - a_4_a_rhs
        self.add_constraints(a_ub_4_a, b_ub_4_a)

        # b) Less than max flow supported by the service station
        ids = [st[0] for st in self.r_to_st]
        r_s_max_i = self.parameters.stations.r_s_max[ids]

        a_ub_4_b = a_4_a_lhs
        b_ub_4_b = np.kron(self.k_ones, r_s_max_i.reshape(-1, 1))
        self.add_constraints(a_ub_4_b, b_ub_4_b)

        # c) Less than max flow prescribed by controller
        if self.controlled:
            a_ub_4_c = a_ub_4_b
            b_ub_4_c = self.control.r_s_c.reshape(-1, 1)
            self.add_constraints(a_ub_4_c, b_ub_4_c)

        # 6) Less than priority times supply
        # a) Less than max flow supported due to congestion in next cell
        j_r = self.parameters.stations.j_r

        p_s_i = np.array([self.parameters.stations.j_to_p[j] for j in j_r])
        p_s_i_k = np.kron(self.k_ones, p_s_i.reshape(-1, 1))
        w_i = self.parameters.highway.w[j_r]
        w_i_k = np.kron(self.k_ones, w_i.reshape(-1, 1))
        rho_max_i = self.parameters.highway.rho_max[j_r]
        rho_max_i_k = np.kron(self.k_ones, rho_max_i.reshape(-1, 1))

        a_6_1 = np.zeros((self.n_r, self.m))
        a_6_1[:, self.n_c + 1:] = np.eye(self.n_r)
        a_6_a_lhs = block_diag(*[a_6_1 for _ in range(self.k_l)])

        r_6 = [(self.n_c * k) + j for k in range(self.k_l) for j in j_r]
        a_6_4 = np.multiply(p_s_i_k, w_i_k)
        a_6_a_rhs = np.multiply(a_6_4, self.rho_matrix[r_6, :])

        a_ub_6_a = a_6_a_lhs + a_6_a_rhs
        b_ub_6_a = np.multiply(a_6_4, (rho_max_i_k - np.matmul(self.rho_0_matrix[r_6, :], self.x_0)))  # Initial condition update
        self.add_constraints(a_ub_6_a, b_ub_6_a)

        # b) Less than max flow supported by the next cell
        q_max_i = self.parameters.highway.q_max[j_r]
        q_max_i_k = np.kron(self.k_ones, q_max_i.reshape(-1, 1))

        a_ub_6_b = a_6_a_lhs
        b_ub_6_b = np.multiply(p_s_i_k, q_max_i_k)
        self.add_constraints(a_ub_6_b, b_ub_6_b)

        # # # #################################
        # #
        # # 5) Less than supply minus mainstream demand
        # a_5_lhs = a_6_a_lhs
        # a_5_rhs_1 = np.multiply(w_i_k, self.rho_matrix[c.tolist(), :])
        # b_ub_5_1 = np.multiply(p_s_i_k, np.multiply(w_i_k, rho_max_i_k))
        # b_ub_5_2 = np.multiply(p_s_i_k, q_max_i_k)
        # beta_ms_j = 1 - self.beta[:, j_r-1]
        # beta_ms_j_k = beta_ms_j.reshape(-1, 1)
        # v_j = self.parameters.highway.v[j_r-1]
        # v_j_k = np.kron(self.k_ones, np.reshape(v_j, (-1, 1)))
        # a_5_1 = np.multiply(beta_ms_j_k, v_j_k)
        # a_5_rhs_3 = np.multiply(a_5_1, self.rho_matrix[(c-1).tolist(), :])  # TODO: Bug here
        # b_ub_5_3 = np.zeros((self.k_int*len(j_r), 1))
        # b_ub_5_4 = np.kron(self.k_ones, self.parameters.highway.q_max[j_r])
        #
        # a_ub_5_a = a_5_lhs + a_5_rhs_1 + a_5_rhs_3
        # b_ub_5_a = b_ub_5_1 - b_ub_5_3
        # self.add_constraints(a_ub_5_a, b_ub_5_a)
        #
        # # a_ub_5_b = a_5_lhs + a_5_rhs_1
        # # b_ub_5_b = b_ub_5_1 - b_ub_5_4
        # # self.add_constraints(a_ub_5_b, b_ub_5_b)
        #
        # a_ub_5_c = a_5_lhs + a_5_rhs_3
        # b_ub_5_c = b_ub_5_2 - b_ub_5_3
        # self.add_constraints(a_ub_5_c, b_ub_5_c)
        #
        # # a_ub_5_d = a_5_lhs + a_5_rhs_3
        # # b_ub_5_d = b_ub_5_2 - b_ub_5_4
        # # self.add_constraints(a_ub_5_d, b_ub_5_d)
        #
        # # 2) Less than supply minus service-station demand
        # a_2_1 = np.zeros((len(j_r), self.m))
        #
        # for i, j in enumerate(j_r):
        #     a_2_1[i, j] += 1
        #
        # a_2_lhs = block_diag(*[a_2_1 for _ in range(self.k_int)])
        #
        # a_2_rhs_1 = a_5_rhs_1
        # b_ub_2_1 = b_ub_5_1
        # b_ub_2_2 = b_ub_5_2
        #
        # a_2_rhs_3 = np.zeros((len(j_r)*self.k_int, self.num_variables))
        #
        # for j in j_r:
        #     r = self.j_to_r[j]
        #     r_1 = [k*len(j_r) + r[0] for k in range(self.k_int)]
        #     r_2 = [k*self.n_r + r[0] for k in range(self.k_int)]
        #     a_2_rhs_3[r_1, :] = a_4_a_rhs[r_2, :]
        #
        # b_ub_2_3 = np.kron(self.k_ones, np.zeros((len(j_r), 1)))
        # r = [self.j_to_r[j] for j in j_r]
        # ids = [self.r_to_st[r_] for r_ in r[0]]
        # r_s_max_j = self.parameters.stations.r_s_max[tuple(ids)]
        # b_ub_2_4 = np.kron(self.k_ones, r_s_max_j)
        #
        # a_ub_2_a = a_2_lhs + a_2_rhs_1 + a_2_rhs_3
        # b_ub_2_a = b_ub_2_1 - b_ub_2_3
        # self.add_constraints(a_ub_2_a, b_ub_2_a)
        #
        # # a_ub_2_b = a_2_lhs + a_2_rhs_1
        # # b_ub_2_b = b_ub_2_1 - b_ub_2_4
        # # self.add_constraints(a_ub_2_b, b_ub_2_b)
        #
        # a_ub_2_c = a_2_lhs + a_2_rhs_3
        # b_ub_2_c = b_ub_2_2 - b_ub_2_3
        # self.add_constraints(a_ub_2_c, b_ub_2_c)
        #
        # # a_ub_2_d = a_2_lhs + a_2_rhs_3
        # # b_ub_2_d = b_ub_2_2 - b_ub_2_4
        # # self.add_constraints(a_ub_2_d, b_ub_2_d)
        #
        # if self.controlled:
        #     b_ub_2_5 = self.control.r_s_c[:, self.j_to_r[j_r]].reshape(-1, 1)
        #     a_ub_2_e = a_2_lhs + a_2_rhs_1
        #     b_ub_2_e = b_ub_2_1 - b_ub_2_5
        #     a_ub_2_f = a_2_lhs
        #     b_ub_2_f = b_ub_2_2 - b_ub_2_5
        #
        #     self.add_constraints(a_ub_2_e, b_ub_2_e)
        #     self.add_constraints(a_ub_2_f, b_ub_2_f)
        #
        # # # #################################

        # Enforce flow into last cell is equivalent to terminal flow
        a_f = np.zeros(self.m)
        a_f[self.n_c - 1] = 1
        a_f[self.n_c] = -1
        a_eq = block_diag(*[a_f for _ in range(self.k_l)])
        b_eq = np.zeros((self.k_l, 1))
        self.add_constraints(a_eq, b_eq, mosek.boundkey.fx)

    # TODO: Test this function!!
    def define_cost_function(self):
        """
        TODO: Add description / comments
        """

        # Define coefficients from Total Time Spent (TTS: Minimize)
        # 1) Minimize highway cell density
        l_i_k = np.reshape(np.kron(self.k_1_ones, self.parameters.highway.l), (-1, 1))
        tts_1 = np.sum(np.multiply(l_i_k, self.rho_matrix), axis=0)

        # 2) Minimize service station exit queue length
        tts_2 = np.sum(self.e_matrix, axis=0)
        # tts = self.dt * (tts_1 + tts_2)

        # 3) Maximize service station use
        tts_3 = 1e9 * np.sum(self.l_matrix, axis=0)
        tts = self.dt * (tts_1 + tts_2 - tts_3)

        # Define coefficients from Total Travel Distance (TTD: Maximize)
        # Maximize flows
        l_1 = np.concatenate((self.parameters.highway.l, np.array([0]), self.parameters.stations.l_r))
        l_2 = np.ones(self.k_l)
        ttd = self.dt * np.kron(l_2, l_1)

        self.c = tts - (self.eta * ttd)

        # Set cost coefficients corresponding to maximize initial flow
        self.c[::self.m] = -10

        # for id in self.parameters.stations.id:
        #     i = self.parameters.stations.i[id]
        #     self.c[i+1::self.m] = -10

    def update_variables(self):
        """
        TODO: Add comments
        """

        self.u_opt_k = np.reshape(self.u_opt, (self.k_l, self.m))

        x_rho = np.matmul(self.rho_0_matrix, self.x_0) + np.matmul(self.rho_matrix, self.u_opt)
        x_l = np.matmul(self.l_0_matrix, self.x_0) + np.matmul(self.l_matrix, self.u_opt) + np.matmul(self.l_s_0_matrix, self.phi_s_0)
        x_e = np.matmul(self.e_0_matrix, self.x_0) + np.matmul(self.e_matrix, self.u_opt) + np.matmul(self.e_s_0_matrix, self.phi_s_0)

        self.x_rho = x_rho.reshape(self.k_l + 1, self.n_c)
        self.x_l = x_l.reshape(self.k_l + 1, self.n_s)
        self.x_e = x_e.reshape(self.k_l + 1, self.n_s)

        self.u_phi = self.u_opt_k[:, :self.n_c + 1]
        self.u_rs = self.u_opt_k[:, self.n_c + 1:]

    def get_control(self) -> TrafficControlInput:
        """
        TODO: !!!
        """
        raise NotImplementedError

    def streamprinter(self, text):
        """
        Define a stream printer to grab output from MOSEK
        """

        sys.stdout.write(text)
        sys.stdout.flush()

    def solve_init(self, x_0: np.ndarray, phi_0: np.ndarray, phi_s_0: np.ndarray, k_0: int, k_l: int):
        """
        Dimension of State: k_0 ~ k_0 + k_window        (Dimension = k_window + 1)
        Dimension of Input: k_0 ~ k_0 + k_window - 1    (Dimension = k_window)
        """

        self.x_0 = x_0
        self.phi_0 = phi_0
        self.phi_s_0 = phi_s_0

        self.k_0 = k_0
        self.k_l = k_l
        k_l_1 = self.k_l + 1
        self.k_ones = np.ones((self.k_l, 1))
        self.k_1_ones = np.ones((k_l_1, 1))

        self.num_variables = self.k_l * self.m
        self.beta = np.zeros((self.k_l, self.n_c))

        self.rho_0_matrix = np.zeros((k_l_1 * self.n_c, self.n))
        self.l_0_matrix = np.zeros((k_l_1 * self.n_s, self.n))
        self.l_s_0_matrix = np.zeros((k_l_1 * self.n_s, len(self.phi_s_0)))
        self.e_0_matrix = np.zeros((k_l_1 * self.n_s, self.n))
        self.e_s_0_matrix = np.zeros((k_l_1 * self.n_s, len(self.phi_s_0)))

        self.rho_matrix = np.zeros((k_l_1 * self.n_c, self.num_variables))
        self.l_matrix = np.zeros((k_l_1 * self.n_s, self.num_variables))
        self.e_matrix = np.zeros((k_l_1 * self.n_s, self.num_variables))

        self.u_opt = np.zeros(self.num_variables)
        self.u_opt_k = np.zeros((self.k_l, self.m))

        self.x_rho = np.zeros((k_l_1, self.n_c))
        self.x_l = np.zeros((k_l_1, self.n_s))
        self.x_e = np.zeros((k_l_1, self.n_s))

        self.u_phi = np.zeros((self.k_l, self.n_c + 1))  # Including terminal flow
        self.u_rs = np.zeros((self.k_l, self.n_r))

        self.c = np.zeros(self.num_variables)

        self.compute_total_beta()
        self.compute_state_matrices()
        self.define_state_constraints()
        self.define_flow_constraints()
        self.formulate_constraints()
        self.define_cost_function()

    def solve(self, print_sol=False):
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

                # Set the bounds on flows
                # blx[j] <= x_j <= bux[j]
                (k, r) = divmod(j, self.m)

                # Initial flow condition
                if r == 0:
                    task.putvarbound(j,  # NOTE: Changed from equality constraint
                                     mosek.boundkey.ra,  # TODO: ...
                                     0,  # TODO: ...
                                     self.phi_0[k])  # TODO: ...
                # Positive flows
                else:
                    task.putvarbound(j,  # TODO: ...
                                     mosek.boundkey.lo,  # TODO: ...
                                     0,  # TODO: ...
                                     +inf)  # TODO: ...

                # Input column j of A (constraints)
                task.putacol(j,  # Variable (column) index.
                             self.asub[j],  # Row index of non-zeros in column j.
                             self.aval[j])  # Non-zero Values of column j.

            # Set the bounds on constraints.
            # b_lb[i] <= constraint_i <= b_ub[i]
            for i in range(self.num_constraints):
                task.putconbound(i,  # TODO: ...
                                 self.bkc[i],  # TODO: ...
                                 self.b_lb[i][0],  # TODO: ...
                                 self.b_ub[i][0])  # TODO: ...

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
                self.u_opt = np.reshape(task.getxx(mosek.soltype.bas), (-1, 1))
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
