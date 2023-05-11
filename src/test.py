import os
import unittest
import numpy as np

from data import read_parameters, read_phi
from optimizer import MyTrafficOptimizer

# TODO: Outdated
class TestReadParameters(unittest.TestCase):
    @classmethod
    def setUp(cls) -> None:
        """Read Parameters and Data"""
        cwd = os.getcwd()
        path = os.path.split(cwd)
        hi_loc = path[0] + '/data/test/cells_1.csv'
        onr_loc = path[0] + '/data/test/onramps_1.csv'
        offr_loc = path[0] + '/data/test/offramps_1.csv'
        stations_loc = path[0] + '/data/test/stations_1.csv'
        phi_loc = path[0] + '/data/phi_1_24h_realsmooth.csv'
        phi_onr_loc = path[0] + '/data/test/onramps_signal_1.csv'

        cls.parameters = read_parameters(hi_loc, onr_loc, offr_loc, stations_loc, phi_onr_loc)
        cls.data = read_phi(phi_loc)

    def test_parameters(self):
        self.assertEqual(len(self.parameters.highway), 4, "incorrect number of cells")

    def test_phi_data(self):
        self.assertEqual(len(self.data), 8640, "incorrect length")

    def test_unique_flows(self):
        correct_flows = np.array([2, 3])
        self.assertTrue(np.allclose(self.parameters.stations.j_r, correct_flows), "incorrect flows")

    def test_priority_update(self):
        correct_priorities = np.array([1, 1, 0.93, 0.83])
        self.assertTrue(np.allclose(self.parameters.highway.p_ms, correct_priorities), "incorrect priorities")

    def test_station_lengths(self):
        correct_lengths = np.array([0.2, 0.4])
        self.assertTrue(np.allclose(self.parameters.stations.l_r, correct_lengths), "incorrect lengths")

    def test_j_to_p(self):
        # j_to_p = {2: 0.93, 3: 0.83}
        # self.assertTrue(self.parameters.stations.j_to_p == j_to_p)
        raise NotImplementedError


class TestOptimizer(unittest.TestCase):
    @classmethod
    def setUp(cls) -> None:
        ####################
        # Input Parameters #
        ####################

        cls.dt = 10                 # Sampling rate [sec]
        cls.window_length = 4       # Length of one interation [hrs]
        cls.eta = 1                 # Optimization Parameter: Tradeoff Between TTS, TTD

        start = 0
        congestion_start = 5
        congestion_end = congestion_start + cls.window_length
        end = 4000

        ############################
        # Read Parameters and Data #
        ############################

        cwd = os.getcwd()
        path = os.path.split(cwd)
        hi_loc = path[0] + '/data/test/cells_1.csv'
        onr_loc = path[0] + '/data/test/onramps_1.csv'
        offr_loc = path[0] + '/data/test/offramps_1.csv'
        stations_loc = path[0] + '/data/test/stations_none.csv'
        phi_loc = path[0] + '/data/test/phi_1.csv'
        phi_onr_loc = path[0] + '/data/test/onramps_signal_1.csv'

        cls.parameters = read_parameters(hi_loc, onr_loc, offr_loc, stations_loc, phi_onr_loc)
        cls.parameters.update_dt(dt=cls.dt/3600)

        cls.opt = MyTrafficOptimizer(parameters=cls.parameters, eta=cls.eta)



        cls.data = read_phi(phi_loc)
        # cls.k_int = 5
        # cls.x0 = np.array([80, 20, 30, 40, 0, 0, 0, 0, 0, 0])
        cls.x0 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


        cls.dt = 0.002777778
        cls.a_i = np.divide(cls.dt, cls.parameters.highway.l)
        cls.beta_i = np.array([0.25, 0.40, 0, 0])
        cls.beta_i_k = np.kron(np.ones((cls.k_int, 1)), cls.beta_i)
        cls.beta_s = np.array([0.05, 0.20, 0.20])
        cls.beta_s_k = np.kron(np.ones((cls.k_int, 1)), cls.beta_s)
        cls.b_i = np.kron(np.ones((cls.k_int, 1)), np.divide(cls.a_i, 1-cls.beta_i))
        cls.b_i_s = np.divide(cls.beta_s_k, 1-cls.beta_i_k[:, [0, 1, 1]]) * cls.dt

        cls.build_state_matrices()

    @classmethod
    def build_state_matrices(cls):
        rho_i_1 = np.array([[cls.a_i[0], -cls.b_i[0, 0], 0, 0, 0, 0],
                            [0, cls.a_i[1], -cls.b_i[0, 1], 0, 0, 0],
                            [0, 0, cls.a_i[2], -cls.b_i[0, 2], cls.a_i[2], 0],
                            [0, 0, 0, cls.a_i[3], 0, cls.a_i[3]]])

        rho_i_2 = np.array([[cls.a_i[0], -cls.b_i[1, 0], 0, 0, 0, 0],
                            [0, cls.a_i[1], -cls.b_i[1, 1], 0, 0, 0],
                            [0, 0, cls.a_i[2], -cls.b_i[1, 2], cls.a_i[2], 0],
                            [0, 0, 0, cls.a_i[3], 0, cls.a_i[3]]])

        rho_i_3 = np.array([[cls.a_i[0], -cls.b_i[2, 0], 0, 0, 0, 0],
                            [0, cls.a_i[1], -cls.b_i[2, 1], 0, 0, 0],
                            [0, 0, cls.a_i[2], -cls.b_i[2, 2], cls.a_i[2], 0],
                            [0, 0, 0, cls.a_i[3], 0, cls.a_i[3]]])

        rho_i_4 = np.array([[cls.a_i[0], -cls.b_i[3, 0], 0, 0, 0, 0],
                            [0, cls.a_i[1], -cls.b_i[3, 1], 0, 0, 0],
                            [0, 0, cls.a_i[2], -cls.b_i[3, 2], cls.a_i[2], 0],
                            [0, 0, 0, cls.a_i[3], 0, cls.a_i[3]]])

        z = np.zeros((4, 6))

        cls.rho_matrix = np.block([[z, z, z, z, z],
                                   [rho_i_1, z, z, z, z],
                                   [rho_i_1, rho_i_2, z, z, z],
                                   [rho_i_1, rho_i_2, rho_i_3, z, z],
                                   [rho_i_1, rho_i_2, rho_i_3, rho_i_4, z]])

        l_i_1 = np.array([[cls.b_i_s[0, 0], 0, 0, 0, -cls.dt, 0],
                          [0, cls.b_i_s[0, 1], 0, 0, 0, -cls.dt],
                          [0, cls.b_i_s[0, 2], 0, 0, 0, -cls.dt]])

        l_i_2 = np.array([[cls.b_i_s[1, 0], 0, 0, 0, -cls.dt, 0],
                          [0, cls.b_i_s[1, 1], 0, 0, 0, -cls.dt],
                          [0, cls.b_i_s[1, 2], 0, 0, 0, -cls.dt]])

        l_i_3 = np.array([[cls.b_i_s[2, 0], 0, 0, 0, -cls.dt, 0],
                          [0, cls.b_i_s[2, 1], 0, 0, 0, -cls.dt],
                          [0, cls.b_i_s[2, 2], 0, 0, 0, -cls.dt]])

        l_i_4 = np.array([[cls.b_i_s[3, 0], 0, 0, 0, -cls.dt, 0],
                          [0, cls.b_i_s[3, 1], 0, 0, 0, -cls.dt],
                          [0, cls.b_i_s[3, 2], 0, 0, 0, -cls.dt]])

        z = np.zeros((3, 6))

        cls.l_matrix = np.block([[z, z, z, z, z],
                                 [l_i_1, z, z, z, z],
                                 [l_i_1, l_i_2, z, z, z],
                                 [l_i_1, l_i_2, l_i_3, z, z],
                                 [l_i_1, l_i_2, l_i_3, l_i_4, z]])

        e_i_1 = np.array([[0, 0, 0, 0, -cls.dt, 0],
                          [0, 0, 0, 0, 0, -cls.dt],
                          [0, 0, 0, 0, 0, -cls.dt]])

        cls.e_matrix = np.block([[z, z, z, z, z],
                             [e_i_1, z, z, z, z],
                             [e_i_1, e_i_1, z, z, z],
                             [e_i_1, e_i_1, e_i_1, z, z],
                             [e_i_1, e_i_1, e_i_1, e_i_1, z]])

        cls.e_matrix[6, :] += cls.l_matrix[3, :]
        cls.e_matrix[7, :] += cls.l_matrix[4, :]

        cls.e_matrix[9, :] += cls.l_matrix[3, :] + cls.l_matrix[6, :]
        cls.e_matrix[10, :] += cls.l_matrix[4, :] + cls.l_matrix[7, :]
        cls.e_matrix[11, :] += cls.l_matrix[5, :]

        cls.e_matrix[12, :] += cls.l_matrix[3, :] + cls.l_matrix[6, :] + cls.l_matrix[9, :]
        cls.e_matrix[13, :] += cls.l_matrix[4, :] + cls.l_matrix[7, :] + cls.l_matrix[10, :]
        cls.e_matrix[14, :] += cls.l_matrix[5, :] + cls.l_matrix[8, :]

        return

    # def test_st_to_r(self):
    #     correct_map = np.array([0, 1, 1])
    #     self.assertTrue(np.all(self.opt.st_to_r == correct_map), "incorrect map")
    #
    # def test_r_to_st(self):
    #     correct_map = [[0], [1, 2]]
    #     self.assertEqual(self.opt.r_to_st, correct_map, "incorrect map")

    def test_station_maps(self):
        raise NotImplementedError

    def test_beta(self):
        self.assertTrue(np.all(self.opt.beta == self.beta_i_k), "incorrect beta")

    def test_beta_s(self):
        self.assertTrue(np.all(self.opt.beta_s == self.beta_s_k), "incorrect beta")

    def test_rho_matrix(self):
        self.assertTrue(np.allclose(self.opt.rho_matrix, self.rho_matrix), "incorrect matrix")

    def test_l_matrix(self):
        self.assertTrue(np.allclose(self.opt.l_matrix, self.l_matrix), "incorrect matrix")

    def test_e_matrix(self):
        self.assertTrue(np.allclose(self.opt.e_matrix, self.e_matrix), "incorrect matrix")

    def test_constraints(self):
        raise NotImplementedError

    def test_formulate_constraints(self):
        raise NotImplementedError

    def test_cost_function(self):
        raise NotImplementedError

    def test_solve(self):
        self.opt.solve(print_sol=True)