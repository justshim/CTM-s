from dataclasses import dataclass
from typing import List, TypeVar

import numpy as np
from matplotlib import pyplot as plt

from model.supervisor import Stretch


P = TypeVar('P', bound='TrafficPerformance')


@dataclass
class TrafficPerformance:
    """
    Dataclass to define, compute performance metrics
    """

    cong_sev: float         # TODO: Congestion Severity
    e_time: float           # Minutes in the Day with Queue [min]
    e_wait: float           # Average Queue Wait Time [min]

    def as_array(self) -> np.ndarray:
        return np.array([self.cong_sev, self.e_time, self.e_wait])

    @staticmethod
    def from_array(p: np.ndarray) -> P:
        return TrafficPerformance(p[0], p[1], p[2])


class TrafficEvaluator:
    """
    Wrapper class to compute TrafficPerformance
    """

    def evaluate(self, s: Stretch, print_perf: bool) -> TrafficPerformance:
        """
        Return TrafficPerformance object
        """

        congestion_severity = self.compute_congestion_severity(s)
        e_time = self.compute_e_time(s)
        e_wait = self.compute_e_wait(s)

        perf = TrafficPerformance(congestion_severity, e_time, e_wait)

        if print_perf:
            print(perf)

        return perf

    @classmethod
    def compute_congestion_severity(cls, s: Stretch):
        """
        TODO: Review the intuition of this
        """

        # Compute rho_critical
        rho_critical = np.zeros(s.n_cells)

        for i, c_i in enumerate(s.cells):
            if i == 0:
                rho_critical[i] = c_i.rho_max - (c_i.q_max / c_i.w)
                continue

            c_prev = s.cells[i-1]
            beta_prev = s.beta_total[0, i-1]  # Assume constant beta

            rho_critical[i] = max(
                [c_i.rho_max - (c_prev.q_max/c_i.w),
                 c_i.w * c_i.rho_max / ((1 - beta_prev) * c_prev.v_free + c_i.w),
                 c_i.rho_max - (c_i.q_max/c_i.w)]
            )

        # Compute congestion severity
        congestion_severity = 0

        for i, c in enumerate(s.cells):
            rho = np.array(c.rho)
            # print(np.sum(rho >= rho_critical[i]))
            rho[rho < rho_critical[i]] = 0
            rho[rho >= rho_critical[i]] -= rho_critical[i]
            congestion_severity += np.sum(np.power(rho, 2))  # TODO: 2-norm?

        return congestion_severity

    @classmethod
    def compute_e_time(cls, s: Stretch):
        """
        Time with queue at the service station (minutes)
        """
        e_time = 0

        for st in s.stations:
            e_time += np.sum(np.array(st.e) > 0)

        return e_time / 6  # Deciseconds to minutes

    @classmethod
    def compute_e_wait(cls, s: Stretch):
        """
        Average wait time in the queue (minutes)
        """
        e_wait = 0

        for st in s.stations:
            e_wait += np.sum(np.array(st.e)) / np.sum(np.array(st.s_e))

        return e_wait / 6  # Deciseconds to minutes


class TrafficResults:
    """
    Class to keep track of TrafficPerformance results
    """

    n_iter: int
    n_fact: int

    def __init__(self, n_iter: int, n_fact: int):
        self.n_iter = n_iter
        self.n_fact = n_fact

    @staticmethod
    def save(j: int, n: int, s: Stretch, p: TrafficPerformance, results_path: str):
        """
        Save state, input arrays from simulation as .npz files
        Save performance result object as .npy file
        """

        np.savez(f"{results_path}f_{j}_it_{n}_cells.npz",
                 rho=s.y_rho,
                 phi=s.u_phi,
                 dem=s.u_dem,
                 sup=s.u_sup,
                 cong=s.u_cong)

        np.savez(f"{results_path}f_{j}_it_{n}_stations.npz",
                 l=s.y_l,
                 e=s.y_e,
                 r_s=s.u_rs)

        np.save(f"{results_path}f_{j}_it_{n}_perf.npy",
                p.as_array())
