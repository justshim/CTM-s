from dataclasses import dataclass

import numpy as np


@dataclass
class ControlParameters:
    """
    TODO: ...
    """

    discount = 0.05                 # TODO: ...
    epi_st = 0.01                    # TODO: ...
    a_rho = 0.003                   # TODO: ...
    a_queue = 0.0009                # TODO: ...
    a_epi = 0.003                   # TODO: ...
    a_cont = 0.70                   # TODO: ...


class TrafficHistory:
    """
    Class to track and predict incoming flows
    """

    n_iter: int                     # Number of iterations
    n_updates: int                  # Number of updates during congestion window
    i: int                          # Index variable
    phi_0_hist: np.ndarray          # Array of all previously observed phi_0
    s: np.random.Generator          # Numpy Generator for flow predictions

    def __init__(self, n_iter: int, window_length: int, max_updates):
        self.n_iter = n_iter
        self.n_updates = max_updates
        self.i = 0
        self.phi_0_hist = np.zeros((window_length, self.n_updates, n_iter))
        self.s = np.random.default_rng(0)

    def log_flow(self, phi_new: np.ndarray):
        """
        Log observed flow into phi_0_hist[~:,j,k]

        ~: length of the observed flow
        j: corresponding to the time instance k_0 when the flow is observed
        k: kth iteration instance of the flow being observed
        """

        k, j = divmod(self.i, self.n_updates)

        self.phi_0_hist[:, j, k] = phi_new
        self.i += 1

    def get_flow(self, predict=False) -> np.ndarray:
        """
        Return random, convex-combination of previously observed flows into the first cell
        Otherwise, return the first observed instance of the flow
        """

        k, j = divmod(self.i, self.n_updates)

        if predict:
            if k > 1:
                s = self.s.dirichlet((k-1, 1))
            else:
                s = 1.00

            return np.sum(np.multiply(s, self.phi_0_hist[:, j, :k]), axis=1)
        else:
            return self.phi_0_hist[:, j, 0]
