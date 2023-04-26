from typing import List

import numpy as np
from numpy.random import Generator


class ControlParameters:
    Q: np.ndarray
    R: np.ndarray

    # Q: np.ndarray
    # R: np.ndarray
    # G: np.ndarray
    # w: np.ndarray
    # d: np.ndarray
    # M: np.ndarray
    #
    # H: np.ndarray
    # f: np.ndarray

    # eta: float

    # def __init__(self):

class TrafficHistory:
    """
    TODO: ...
    """

    n_iter: int                     # Number of Iterations
    n_updates: int                  # Number of updates during congestion window
    i: int                          # Index variable
    phi_0_hist: np.ndarray          # TODO: Array of all previously observed phi_0
    s: Generator                    # TODO: ...

    def __init__(self, n_iter: int, window_length: int, max_updates):
        self.n_iter = n_iter
        self.n_updates = max_updates
        self.i = 0
        self.phi_0_hist = np.zeros((window_length, self.n_updates, n_iter))
        self.s = np.random.default_rng(0)

    def log_flow(self, phi_new: np.ndarray):
        """
        TODO: ...
        """

        k, j = divmod(self.i, self.n_updates)

        self.phi_0_hist[:, j, k] = phi_new
        self.i += 1

    def estimate_flow(self) -> np.ndarray:
        """
        Return random, convex-combination of previously observed flows into the first cell
        """

        k, j = divmod(self.i, self.n_updates)

        if k > 1:
            s = self.s.dirichlet((k-1, 1))
        else:
            s = 1.00

        return np.sum(np.multiply(s, self.phi_0_hist[:, j, :k]), axis=1)
