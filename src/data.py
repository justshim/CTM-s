import numpy as np
import pandas as pd


class TrafficFlowGenerator:
    """
    TODO: ...
    """
    s: np.random.Generator              # TODO: ...
    phi: np.ndarray                     # TODO: ...

    mu: int
    std: int
    k: int

    n_spikes: int

    def __init__(self, loc: str, ds=1):
        self.s = np.random.default_rng(0)
        self.phi = pd.read_csv(loc, header=None).to_numpy().flatten()

        self.mu = 0
        self.std = 1
        self.k = 50

        self.n_spikes = 10

    def get_flow(self, t_0: int, t_f: int, perturb=False) -> np.ndarray:
        """
        TODO: ...
        """

        if perturb:
            return self.perturb_flow()
        else:
            return self.phi[t_0:t_f]

    def perturb_flow(self) -> np.ndarray:
        """
        TODO: ...
        """

        phi_noisy = np.copy(self.phi)

        # Introduce random Gaussian Noise
        phi_noisy += self.k * self.s.normal(self.mu, self.std, self.phi.shape)

        # Introduce random spikes into the flow
        # TODO: Basically need to ensure that the spikes last for a certain amount of time and are well-spaced
        #  Also need to determine the magnitude of the peaks...
        #
        # np.random.rand(0, len(phi), num_peaks)

        return phi_noisy

def read_phi(loc: str, t_0: int, t_f: int, ds=1) -> np.ndarray:
    """
    Read in real-flow data from t_0 to t_f
    Flow data is measured every 10 seconds over the course of 24 hours
    ds = down-sampling rate e.g. ds = 6 => flow every minute
    """

    phi = pd.read_csv(loc, header=None).to_numpy().flatten()
    return phi[t_0:t_f:int(ds)]


# TODO: Function to write synthetic/perturbed phi_data
def perturb_flow(phi: np.ndarray) -> np.ndarray:
    """
    TODO: ...
    """

    # Introduce random Gaussian Noise
    # TODO: Seed the random generator?
    # TODO: Gaussian Noise Parameters
    mu = 0
    std = 1
    k = 100
    phi += k * np.random.normal(mu, std, phi.shape)

    # Introduce random spikes into the flow
    # TODO: Basically need to ensure that the spikes last for a certain amount of time and are well-spaced
    #  Also need to determine the magnitude of the peaks...
    #
    num_peaks = 10  # TODO: Hardcoded value here
    # np.random.rand(0, len(phi), num_peaks)

    return phi

# TODO: Function to write synthetic service station beta split ratio that will need to be learned via ILC...



