import sys
from typing import List, Dict

import numpy as np
from scipy.linalg import block_diag
from math import inf
import mosek

from model.parameters import CTMsParameters
from control import ControlParameters
from results import plot_lp


# TODO: Cleanup use of np.reshape() throughout the class
class TrafficOptimizer:
    """
    Class to define optimization problem
    """

    params: CTMsParameters          # CTM-s Parameters Object
    params_c: ControlParameters     # Traffic Control Parameters Object
    x_0: np.ndarray                 # Initial Highway State at time of LP initialization
    phi_0: np.ndarray               # Incoming Flow at First Cell (Estimate)
    s_s_0: np.ndarray               # Service Station in-flow history

    k_0: int                        # Time at LP initialization
    k_l: int                        # Length of LP solving horizon [# of time increments]
    dt: float                       # Time Increment [hrs]

    n_c: int                        # Number of Cells
    n_on: int                       # Number of On-ramps
    n_off: int                      # Number of Off-ramps
    n_s: int                        # Number of Stations
    n_r: int                        # Number of Station Exit Flow Points
    n: int                          # State Dimension: Cell density, number of users and queue length for each station
    m: int                          # Input Dimension: Flows into each cell, flows exiting stations, terminal flow

    beta: np.ndarray                # Total Split Ratio of each Cell
    beta_s: np.ndarray              # Service Station Split Ratio of each Cell

    rho_0_matrix: np.ndarray        # [Dynamics]: Initial Condition Matrix for Cell Density
    l_0_matrix: np.ndarray          # [Dynamics]: Initial Condition Matrix for Number of Station Users
    l_s_0_matrix: np.ndarray        # [Dynamics]: Service Station Inflow Matrix for Number of Station Users
    e_0_matrix: np.ndarray          # [Dynamics]: Initial Condition Matrix for Station Queue Length
    e_s_0_matrix: np.ndarray        # [Dynamics]: Service Station Inflow Matrix for Queue Length
    rho_matrix: np.ndarray          # [Dynamics]: Input to State Matrix (Flows to Cell Density)
    l_matrix: np.ndarray            # [Dynamics]: Input to State Matrix (Flows to Number of Station Users)
    e_matrix: np.ndarray            # [Dynamics]: Input to State Matrix (Flows to Station Queue Length)

    num_constraints: int            # [Mosek]: Number of Constraints in LP
    num_variables: int              # [Mosek]: Number of Optimization Variables in LP
    a_ub: List                      # [Mosek]: Constraint Matrix
    asub: List                      # [Mosek]: Non-zero indices of columns in constraint matrix
    aval: List                      # [Mosek]: Non-zero values of columns in constraint matrix
    b_lb: List                      # [Mosek]: Lower Constraint
    b_ub: List                      # [Mosek]: Upper Constraint
    bkc: List                       # [Mosek]: Constraint Bound Keys
    c: np.ndarray                   # [Mosek]: Cost Function Matrix

    c_opt: float                    # [Solution]: Optimal Cost

    u_opt: np.ndarray               # [Solution]: Optimal Input
    u_opt_k: np.ndarray             # [Solution]: Optimal Input, each time instant k
    u_phi: np.ndarray               # [Solution]: Traffic Flow Between Cells (Input)
    u_rs_c: np.ndarray              # [Solution]: Traffic Flow from Station (Input/Control)

    y_rho: np.ndarray               # [Solution]: Cell Density (Output)
    y_l: np.ndarray                 # [Solution]: Number of People at Station (Output)
    y_e: np.ndarray                 # [Solution]: Station Queue Length (Output)

    st_to_r: np.ndarray             # [Helper]: Station ID to Input index (after cell flows)
    r_to_st: List                   # [Helper]: Input index (after cell flows) to Station ID
    j_to_r: Dict                    # [Helper]: Station outlet cell to index of input (after cell flows)
    k_ones: np.ndarray              # [Helper]: Vector of 1s with length k_l
    k_1_ones: np.ndarray            # [Helper]: Vector of 1s with length k_l+1
    discount: np.ndarray            # [Helper]: Exponential decay profile for cost function

    def __init__(self, params: CTMsParameters, params_c: ControlParameters):
        self.params = params
        self.params_c = params_c
        self.beta_s = self.params.stations.beta_s
        self.dt = params.highway.dt

        self.n_c = len(self.params.highway)
        self.n_on = len(self.params.onramps)
        self.n_off = len(self.params.offramps)
        self.n_s = len(self.params.stations)
        self.n_r = len(self.params.stations.j_r)
        self.n = self.n_c + (2 * self.n_s)
        self.m = self.n_c + self.n_r + self.n_on + 1

        self.build_station_outflow_mapping()

    def build_station_outflow_mapping(self):
        """
        Build maps from/to station index,station outlet cell, and input variable index
        """

        self.st_to_r = np.zeros_like(self.params.stations.j)
        self.r_to_st = [[] for _ in self.params.stations.j_r]
        self.j_to_r = {}

        for i, j in enumerate(self.params.stations.j):
            for r_, j_r in enumerate(self.params.stations.j_r):
                if j == j_r:
                    self.st_to_r[i] = r_
                    self.r_to_st[r_].append(i)
                    self.j_to_r[j] = np.array([r_])
                    break

    def compute_total_beta(self):
        """
        Compute total beta from service-stations and off-ramps
        """

        if self.params.stations.beta_s.ndim == 1:
            self.beta_s = np.kron(self.k_ones, self.params.stations.beta_s)

        for k in range(self.k_l):
            for i_, i in enumerate(self.params.stations.i):
                self.beta[k, i] += self.beta_s[k, i_]

            for i_, i in enumerate(self.params.offramps.i):
                self.beta[k, i] += self.params.offramps.beta_r[i_]

    def compute_state_matrices(self):
        """
        Computation of state dynamic matrices in lifted system representation

        rho_i(k) = (rho_0_matrix * initial state) + (rho_matrix * input)
        l_q(k) = (l_0_matrix * initial state) + (l_s_0_matrix * station inflow history) + (l_matrix * input)
        e_q(k) = (e_0_matrix * initial state) + (e_s_0_matrix * station inflow history) + (e_matrix * input)
        """

        a_i = np.divide(self.dt, self.params.highway.l)
        a_i_k = np.kron(self.k_ones, a_i)
        b_i = np.divide(a_i_k, (1 - self.beta))
        b_i_s = self.beta_s / (1 - self.beta[:, self.params.stations.i])

        #############################
        # Cell Density (Rho Matrix) #
        #############################

        # Initial condition matrix
        self.rho_0_matrix[:, :self.n_c] = np.kron(self.k_1_ones, np.eye(self.n_c))

        # Input to state matrix
        rho_1 = np.diag(a_i)
        rho_4 = np.zeros((self.n_c, self.n_r))

        for k in range(1, self.k_l + 1):
            rho_2 = np.diag(b_i[k - 1, :])
            rho_3 = np.zeros((self.n_c, self.n_c + 1))
            rho_3[:, :self.n_c] += rho_1
            rho_3[:, 1:self.n_c + 1] -= rho_2
            rho_5 = np.concatenate((rho_3, rho_4), axis=1)

            for i_, j_r in enumerate(self.params.stations.j_r):
                rho_5[j_r, self.n_c + 1 + self.st_to_r[i_]] += a_i[j_r]

            rho_6 = np.ones((self.k_l - k + 1, 1))
            self.rho_matrix[k * self.n_c:, (k - 1) * self.m:k * self.m] = np.kron(rho_6, rho_5)

        ###############################################
        # Service Station Number of People (L Matrix) #
        ###############################################

        # Initial condition matrix
        l_0 = np.zeros((self.n_s, self.n))
        for i in self.params.stations.id:
            l_0[i, self.n_c + (2 * i)] = 1

        self.l_0_matrix = np.kron(self.k_1_ones, l_0)

        # Input to state matrix
        # Inflow term
        for k in range(1, self.k_l + 1):
            l_1 = np.zeros((self.n_s, self.m))

            for i_, i in enumerate(self.params.stations.i):
                l_1[i_, i + 1] += self.dt * b_i_s[k - 1, i_]

            l_2 = np.ones((self.k_l - k + 1, 1))
            l_3 = np.kron(l_2, l_1)
            self.l_matrix[k * self.n_s:, (k - 1) * self.m:k * self.m] += l_3

        # Outflow term
        for k in range(1, self.k_l + 1):
            for i_, dq in enumerate(self.params.stations.delta):
                i = self.params.stations.i[i_]
                r_1 = (k * self.n_s) + i_

                k_dq_1 = self.k_0 - round(dq)
                k_dq_2 = self.k_0 - round(dq) + k

                if k_dq_2 > 0:
                    for k_dq in range(max(k_dq_1, 0), min(k_dq_2, self.k_0)):
                        self.l_s_0_matrix[r_1, k_dq] -= self.dt  # b_i_s included in self.s_s_0  # TODO

                    if k_dq_2 > self.k_0:
                        for k_dq in range(self.k_0, k_dq_2):
                            k_dq_3 = k_dq - self.k_0
                            c_1 = (k_dq_3 * self.m) + (i + 1)
                            self.l_matrix[r_1, c_1] -= (self.dt * b_i_s[k_dq_3, i_])  # TODO

        ###########################################
        # Service Station Queue Length (E Matrix) #
        ###########################################

        # Initial condition matrix
        e_0 = np.zeros((self.n_s, self.n))
        for i in self.params.stations.id:
            e_0[i, self.n_c + 1 + (2 * i)] += 1

        self.e_0_matrix = np.kron(self.k_1_ones, e_0)

        # Input to state matrix
        # Outflow term
        for k in range(1, self.k_l + 1):
            e_1 = np.zeros((self.n_s, self.m))
            for i in self.params.stations.id:
                e_1[i, self.n_c + 1 + self.st_to_r[i]] -= self.dt

            e_2 = np.ones((self.k_l - k + 1, 1))
            e_3 = np.kron(e_2, e_1)
            self.e_matrix[k * self.n_s:, (k - 1) * self.m:k * self.m] += e_3

        # Inflow term
        for k in range(1, self.k_l + 1):
            for i_, dq in enumerate(self.params.stations.delta):
                i = self.params.stations.i[i_]
                r_1 = (k * self.n_s) + i_

                k_dq_1 = self.k_0 - round(dq)
                k_dq_2 = self.k_0 - round(dq) + k

                if k_dq_2 > 0:
                    for k_dq in range(max(k_dq_1, 0), min(k_dq_2, self.k_0)):
                        self.e_s_0_matrix[r_1, k_dq] += self.dt  # b_i_s included in self.s_s_0
                    if k_dq_2 > self.k_0:
                        for k_dq in range(self.k_0, k_dq_2):
                            k_dq_3 = k_dq - self.k_0
                            c_1 = (k_dq_3 * self.m) + (i + 1)
                            self.e_matrix[r_1, c_1] += (self.dt * b_i_s[k_dq_3, i_])  # TODO

    def add_constraints(self, a: np.ndarray, b: np.ndarray, bkc=mosek.boundkey.up):
        """
        Helper function to build constraints for Mosek solver
        """

        assert (a.shape[0] == b.shape[0])

        self.a_ub += a.tolist()

        if bkc == mosek.boundkey.ra:
            self.b_lb += b[:, 0].reshape(-1, 1).tolist()
            self.b_ub += b[:, 1].reshape(-1, 1).tolist()
        else:
            self.b_lb += b.tolist()
            self.b_ub += b.tolist()

        self.bkc += [bkc for _ in range(len(b))]

    def formulate_constraints(self):
        """
        Helper function to formulate constraint matrix for Mosek solver
        """

        self.num_constraints = len(self.a_ub)
        a_ub = np.array(self.a_ub)

        for column in a_ub.T:
            ind = np.nonzero(column)[0].tolist()
            self.asub += [ind]
            self.aval += [column[ind].tolist()]

    def define_state_constraints(self):
        """
        Define constraints on cell density, number of station users, station queue length
        """

        ################################
        # Max Cell Density Constraints #
        ################################

        rho_max = self.params.highway.rho_max.reshape((-1, 1))
        b_rho_0 = np.zeros(((self.k_l + 1) * self.n_c, 1)) - np.matmul(self.rho_0_matrix, self.x_0)
        b_rho_max = np.kron(self.k_1_ones, rho_max) - np.matmul(self.rho_0_matrix, self.x_0)
        b_rho = np.concatenate((b_rho_0, b_rho_max), axis=1)
        self.add_constraints(self.rho_matrix, b_rho, mosek.boundkey.ra)

        ############################
        # Queue Length Constraint #
        ############################

        e_max = self.params.stations.e_max.reshape(-1, 1)
        b_e_0 = np.zeros(((self.k_l + 1) * self.n_s, 1)) - np.matmul(self.e_0_matrix, self.x_0) - np.matmul(self.e_s_0_matrix, self.s_s_0)
        b_e_max = np.kron(self.k_1_ones, e_max) - np.matmul(self.e_0_matrix, self.x_0) - np.matmul(self.e_s_0_matrix, self.s_s_0)
        b_e = np.concatenate((b_e_0, b_e_max), axis=1)
        self.add_constraints(self.e_matrix, b_e, mosek.boundkey.ra)

        ############################################
        # Max Service Station Capacity Constraints #
        ############################################
        # TODO: Add max service station capacity constraint...
        b_l_0 = np.zeros(((self.k_l + 1) * self.n_s, 1)) - np.matmul(self.l_0_matrix, self.x_0) - np.matmul(self.l_s_0_matrix, self.s_s_0)
        self.add_constraints(self.l_matrix, b_l_0, mosek.boundkey.lo)

    def define_flow_constraints(self):
        """
        Define constraints on mainstream flows and station out-flows
        """

        ###############################
        # Mainstream Flow Constraints #
        ###############################

        # 1) Less than previous demand
        # a) Less than max flow supported from free-flow in previous cell
        # Note: But no constraint on phi_0, phi_final
        beta_ms_i = 1 - self.beta[:, :-1]
        beta_ms_i_k = beta_ms_i.reshape(-1, 1)
        v_i = self.params.highway.v[:-1]
        v_i_k = np.kron(self.k_ones, v_i.reshape(-1, 1))

        a_1_1 = np.zeros((self.n_c - 1, self.m))
        a_1_1[:, 1:self.n_c] = np.eye(self.n_c - 1)
        a_1_a_lhs = np.zeros((self.k_l * (self.n_c - 1), self.num_variables))
        a_1_a_lhs[:, :-1] = block_diag(*[a_1_1 for _ in range(self.k_l)])
        a_1_2 = np.multiply(beta_ms_i_k, v_i_k)
        r_1 = [(self.n_c * k) + i for k in range(self.k_l) for i in range(self.n_c - 1)]
        a_1_a_rhs = np.multiply(a_1_2, self.rho_matrix[r_1, :])

        a_ub_1_a = a_1_a_lhs - a_1_a_rhs
        b_ub_1_a = np.multiply(a_1_2, np.matmul(self.rho_0_matrix[r_1, :], self.x_0))

        # b) Less than max flow supported by the previous cell
        # Note: But no constraint on phi_0, phi_final
        q_max_i = self.params.highway.q_max[:-1]

        a_ub_1_b = a_1_a_lhs
        b_ub_1_b = np.kron(self.k_ones, q_max_i.reshape(-1, 1))

        # Add all mainstream flow constraints
        self.add_constraints(a_ub_1_a, b_ub_1_a)
        self.add_constraints(a_ub_1_b, b_ub_1_b)

        ####################################
        # Service Station Flow Constraints #
        ####################################

        # 4) Less than previous demand
        # a) Less than max flow possible from service station
        a_4_1 = np.zeros((self.n_r, self.m))
        a_4_1[:, self.n_c + 1:] = np.eye(self.n_r)
        a_4_a_lhs = np.zeros((self.k_l * self.n_r, self.num_variables))
        a_4_a_lhs[:, :-1] = block_diag(*[a_4_1 for _ in range(self.k_l)])
        a_4_a_rhs = np.zeros_like(a_4_a_lhs)
        b_ub_4_a = np.zeros((self.k_l * self.n_r, 1))

        b_i_s = self.beta_s / (1 - self.beta[:, self.params.stations.i])

        for k in range(self.k_l):
            for r_, stations in enumerate(self.r_to_st):
                r_1 = (k * self.n_r) + r_
                for i_ in stations:
                    r_2 = (k * self.n_s) + i_
                    a_4_a_rhs[r_1, :] += (self.e_matrix[r_2, :] / self.dt)  # TODO
                    b_ub_4_a[r_1] += (np.matmul(self.e_0_matrix[r_2, :], self.x_0) / self.dt)
                    b_ub_4_a[r_1] += (np.matmul(self.e_s_0_matrix[r_2, :], self.s_s_0) / self.dt)  # TODO

                    dq = self.params.stations.delta[i_]
                    k_dq_1 = self.k_0 - round(dq)
                    k_dq_2 = k_dq_1 + k

                    i = self.params.stations.i[i_]
                    if k_dq_2 >= 0:
                        if k_dq_2 < self.k_0:
                            b_ub_4_a[r_1] += self.s_s_0[k_dq_2]
                        else:
                            k_dq_3 = k_dq_2 - self.k_0
                            c_1 = (k_dq_3 * self.m) + i + 1
                            a_4_a_rhs[r_1, c_1] += b_i_s[k_dq_3, i_]

        a_ub_4_a = a_4_a_lhs - a_4_a_rhs

        # b) Less than max flow supported by the service station
        ids = [st[0] for st in self.r_to_st]
        r_s_max_j = self.params.stations.r_s_max[ids]

        a_ub_4_b = a_4_a_lhs
        b_ub_4_b = np.kron(self.k_ones, r_s_max_j.reshape(-1, 1))

        # Add all station outflow constraints
        self.add_constraints(a_ub_4_a, b_ub_4_a)
        self.add_constraints(a_ub_4_b, b_ub_4_b)

        #######################################
        # Total Demand and Supply Constraints #
        #######################################

        # Less than supply minus service-station demand
        # Note: For flows unaffected by station outflow, equivalent to ensuring flow is less than supply
        a_3_1 = np.zeros((self.n_c, self.m))
        a_3_1[:, :self.n_c] = np.eye(self.n_c)
        c_3 = [self.n_c + 1 + self.j_to_r[j_r][0] for j_r in self.params.stations.j_r]
        a_3_1[self.params.stations.j_r, c_3] += 1
        a_diff_lhs = np.zeros((self.k_l * self.n_c, self.num_variables))
        a_diff_lhs[:, :-1] = block_diag(*[a_3_1 for _ in range(self.k_l)])

        w_i = self.params.highway.w
        w_i_k = np.kron(self.k_ones, w_i.reshape(-1, 1))
        rho_max_i = self.params.highway.rho_max
        rho_max_i_k = np.kron(self.k_ones, rho_max_i.reshape(-1, 1))
        q_max_i = self.params.highway.q_max
        q_max_i_k = np.kron(self.k_ones, q_max_i.reshape(-1, 1))
        r_3 = [(self.n_c * k) + i for k in range(self.k_l) for i in range(self.n_c)]

        a_diff_rhs = np.multiply(w_i_k, self.rho_matrix[r_3, :])

        b_diff_1 = np.multiply(w_i_k, (rho_max_i_k - np.matmul(self.rho_0_matrix[r_3, :], self.x_0)))
        b_diff_2 = q_max_i_k

        a_diff_1 = a_diff_lhs + a_diff_rhs
        a_diff_2 = a_diff_lhs

        self.add_constraints(a_diff_1, b_diff_1)
        self.add_constraints(a_diff_2, b_diff_2)

        ########################
        # Epigraph Constraints #
        ########################

        a_epigraph = -np.ones(self.num_variables)

        for i in range(self.n_r):
            a_epigraph[self.n_c + 1 + i::self.m] *= self.params_c.epi_st

        a_epigraph[:-1] = np.multiply(a_epigraph[:-1], self.discount)
        a_epigraph = np.atleast_2d(a_epigraph)
        b_epigraph = np.zeros((1, 1))

        self.add_constraints(a_epigraph, b_epigraph)

    def define_cost_function(self):
        """
        TODO: Add description / comments
        """

        # Define control objective
        # Minimize highway cell density [veh/km]
        # Note: Minimizing density rather than number of vehicles
        cost_rho = self.params_c.a_rho * np.sum(self.rho_matrix, axis=0)

        # Minimize service station exit queue length
        cost_queue = self.params_c.a_queue * np.sum(self.e_matrix, axis=0)

        cost = cost_rho + cost_queue  # TODO: Units: [veh hr / km]

        # Exponential discount
        cost[:-1] = np.multiply(cost[:-1], self.discount)

        # Define epigraph of flows
        epi = np.zeros(self.num_variables)
        epi[-1] = 1

        self.c = (self.params_c.a_epi * epi) + (self.params_c.a_cost * cost)

    def update_variables(self):
        """
        Update variables using solution from LP
        """

        self.c_opt = np.matmul(self.c, self.u_opt)[0]
        self.u_opt_k = np.reshape(self.u_opt[:-1, :], (self.k_l, self.m))

        x_rho = np.matmul(self.rho_0_matrix, self.x_0) \
            + np.matmul(self.rho_matrix, self.u_opt)
        self.y_rho = x_rho.reshape(self.k_l + 1, self.n_c)

        x_l = np.matmul(self.l_0_matrix, self.x_0) \
            + np.matmul(self.l_matrix, self.u_opt) \
            + np.matmul(self.l_s_0_matrix, self.s_s_0)
        self.y_l = x_l.reshape(self.k_l + 1, self.n_s)

        x_e = np.matmul(self.e_0_matrix, self.x_0) \
            + np.matmul(self.e_matrix, self.u_opt) \
            + np.matmul(self.e_s_0_matrix, self.s_s_0)
        self.y_e = x_e.reshape(self.k_l + 1, self.n_s)

        self.u_phi = self.u_opt_k[:, :self.n_c + 1]
        self.u_rs_c = self.u_opt_k[:, self.n_c + 1:]

    def solve_init(self, x_0: np.ndarray, phi_0: np.ndarray, s_s_0: np.ndarray, k_0: int, k_l: int):
        """
        Dimension of State: k_0 ~ k_0 + k_l        (Dimension = k_l + 1)
        Dimension of Input: k_0 ~ k_0 + k_l - 1    (Dimension = k_l)
        """

        self.x_0 = x_0
        self.phi_0 = phi_0
        self.s_s_0 = s_s_0

        self.k_0 = k_0
        self.k_l = k_l
        k_l_1 = self.k_l + 1
        self.k_ones = np.ones((self.k_l, 1))
        self.k_1_ones = np.ones((k_l_1, 1))
        self.discount = np.kron(np.exp(-self.params_c.discount * np.arange(0, self.k_l)), np.ones(self.m))

        self.num_variables = (self.k_l * self.m) + 1
        self.beta = np.zeros((self.k_l, self.n_c))

        self.rho_0_matrix = np.zeros((k_l_1 * self.n_c, self.n))
        self.l_0_matrix = np.zeros((k_l_1 * self.n_s, self.n))
        self.l_s_0_matrix = np.zeros((k_l_1 * self.n_s, len(self.s_s_0)))
        self.e_0_matrix = np.zeros((k_l_1 * self.n_s, self.n))
        self.e_s_0_matrix = np.zeros((k_l_1 * self.n_s, len(self.s_s_0)))

        self.rho_matrix = np.zeros((k_l_1 * self.n_c, self.num_variables))
        self.l_matrix = np.zeros((k_l_1 * self.n_s, self.num_variables))
        self.e_matrix = np.zeros((k_l_1 * self.n_s, self.num_variables))

        self.u_opt = np.zeros(self.num_variables)
        self.u_opt_k = np.zeros((self.k_l, self.m))

        self.y_rho = np.zeros((k_l_1, self.n_c))
        self.y_l = np.zeros((k_l_1, self.n_s))
        self.y_e = np.zeros((k_l_1, self.n_s))

        self.u_phi = np.zeros((self.k_l, self.n_c + 1))  # Including terminal flow
        self.u_rs_c = np.zeros((self.k_l, self.n_r))

        self.c = np.zeros(self.num_variables)

        self.a_ub = []
        self.asub = []
        self.aval = []
        self.b_lb = []
        self.b_ub = []
        self.bkc = []

        # TODO: Don't have to recompute all of these every time
        self.compute_total_beta()
        self.compute_state_matrices()
        self.define_state_constraints()
        self.define_flow_constraints()
        self.formulate_constraints()
        self.define_cost_function()

    def solve(self, print_sol=False, plot_sol=False, n_update=0):
        """
        Solve LP for optimal flows using Mosek
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

                if k < self.k_l:
                    # Initial flow condition
                    if r == 0:
                        task.putvarbound(j,                     # Variable (column) index.
                                         mosek.boundkey.ra,     # Lower and Upper Bound
                                         0,                     # Non-negative flow lower bound
                                         self.phi_0[k])         # phi_0 upper bound (from data)
                    # Non-negative flow condition
                    else:
                        task.putvarbound(j,                     # Variable (column) index.
                                         mosek.boundkey.lo,     # Lower bound
                                         0,                     # Non-negative flow lower bound
                                         +inf)
                else:
                    task.putvarbound(j,                         # Variable (column) index.
                                     mosek.boundkey.fr,         # Lower bound
                                     0,
                                     0)

                # Input column j of A (constraints)
                task.putacol(j,                             # Variable (column) index.
                             self.asub[j],                  # Row index of non-zeros in column j.
                             self.aval[j])                  # Non-zero Values of column j.

            # Set the bounds on constraints.
            # b_lb[i] <= constraint_i <= b_ub[i]
            for i in range(self.num_constraints):
                task.putconbound(i,                         # Constraint (row) index
                                 self.bkc[i],               # Mosek Bound Key Constraint Type
                                 self.b_lb[i][0],           # Lower Bound
                                 self.b_ub[i][0])           # Upper Bound

            # Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.minimize)

            # Solve the problem
            task.optimize()

            # Print a summary containing information
            # about the solution for debugging purposes
            if print_sol:
                task.solutionsummary(mosek.streamtype.msg)

            # Get status information about the solution
            solsta = task.getsolsta(mosek.soltype.bas)

            if solsta == mosek.solsta.optimal:
                self.u_opt = np.reshape(task.getxx(mosek.soltype.bas), (-1, 1))

                self.update_variables()
                # self.evaluate_performance()

                if print_sol:
                    print("Optimal solution: ")
                    for i in range(self.num_variables):
                        print("x[" + str(i) + "]=" + str(self.u_opt[i]))

                if plot_sol:
                    plot_lp(self.u_rs_c, [0], self.k_0, self.k_0 + self.k_l, n_update)

            elif (solsta == mosek.solsta.dual_infeas_cer or
                  solsta == mosek.solsta.prim_infeas_cer):
                print("Primal or dual infeasibility certificate found.\n")
            elif solsta == mosek.solsta.unknown:
                print("Unknown solution status")
            else:
                print("Other solution status")

    @staticmethod
    def streamprinter(text):
        """
        Define a stream printer to grab output from MOSEK
        """

        sys.stdout.write(text)
        sys.stdout.flush()
