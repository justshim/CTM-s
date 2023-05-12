import numpy as np
import pandas as pd


# TODO: Play around with noise parameters, function all works though
class TrafficFlowGenerator:
    """
    """
    s: np.random.Generator              # TODO: ...
    phi: np.ndarray                     # TODO: ...

    t_0: int                            # Congestion Start
    t_f: int                            # Congestion End

    mu: int                             # TODO: ...
    std: int                            # TODO: ...
    k_noise: int                        # TODO: ...

    n_spikes: int                       # Number of Spikes
    t_spike: int                        # Duration of the Spike
    k_spike: int                        # Magnitude of Spike
    spike: np.ndarray                   # Spike Profile
    t_spikes: np.ndarray                # Location of Spikes

    def __init__(self, loc: str, t_0: int, t_f: int, ds=1):
        self.s = np.random.default_rng(0)
        self.phi = pd.read_csv(loc, header=None).to_numpy().flatten()

        self.t_0 = t_0
        self.t_f = t_f

        self.mu = 0
        self.std = 1
        self.k_noise = 50

        self.n_spikes = 5
        self.t_spike = 61
        self.k_spike = 200

        self.spike = np.zeros(self.t_spike)
        t_spike_half = int((self.t_spike + 1) / 2)
        self.spike[:t_spike_half] = np.linspace(0, self.k_spike, t_spike_half)
        self.spike[t_spike_half-1:] = np.linspace(self.k_spike, 0, t_spike_half)

        self.t_peaks = np.random.choice(int((self.t_f - self.t_0) / self.t_spike), self.n_spikes, replace=False)

    def get_flow(self, t_0: int, t_f: int, perturb=False) -> np.ndarray:
        """
        Return either nominal or perturbed flow incoming into first cell
        """

        if perturb:
            phi_perturbed = self.perturb_flow()
            return phi_perturbed[t_0:t_f]
        else:
            return self.phi[t_0:t_f]

    def perturb_flow(self) -> np.ndarray:
        """
        Perturb flow profile with

        1) Random Gaussian Noise
        2) Flow spikes at random, pre-determined time instances
        """

        phi_noisy = np.copy(self.phi)

        # Introduce random Gaussian Noise
        phi_noisy += self.k_noise * self.s.normal(self.mu, self.std, self.phi.shape)

        # Introduce random spikes into the flow (same for all flows)
        for t in self.t_peaks:
            t_0 = self.t_0 + (self.t_spike * t)
            t_f = t_0 + self.t_spike
            phi_noisy[t_0:t_f] += self.spike

        return phi_noisy
