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
                c_i.rho_max - (c_prev.q_max/c_i.w),
                c_i.w * c_i.rho_max / ((1 - beta_prev) * c_prev.v_free + c_i.w),
                c_i.rho_max - (c_i.q_max/c_i.w)
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

    def save(self, j: int, n: int, s: Stretch, p: TrafficPerformance, results_path: str):
        """
        Save state, input arrays from simulation as .npz files
        Save performance result object as .npy file
        """

        for i, c in enumerate(s.cells):
            np.savez(f"{results_path}f_{j}_it_{n}_cell_{i}.npz",
                     rho=np.array(c.rho),
                     phi=np.array(c.phi))

        for i, st in enumerate(s.stations):
            np.savez(f"{results_path}f_{j}_it_{n}_stat_{i}.npz",
                     l=np.array(st.l),
                     e=np.array(st.e),
                     r_s=np.array(st.r_s))

        np.save(f"{results_path}f_{j}_it_{n}_perf.npy",
                p.as_array())


# TODO: Get rid of these plotting help functions eventually
def plot_comparison_test(y_0: np.ndarray, y_1: np.ndarray, i: List, k_0: int, k_f: int, loc: str):
    fig = plt.figure()
    t = np.arange(k_0, k_f)
    plt.plot(t, y_0[:, i], '-k')
    plt.plot(t, y_1[:, i], '--r')
    plt.xlabel('Time')
    fig.savefig(loc + 'results_test/cell.png', dpi=300)  # TODO: Naming scheme for files...
    plt.close()


def plot_flow_comparison(y_0: np.ndarray, y_1: np.ndarray, k_0: int, k_f: int, loc: str):
    fig = plt.figure()
    t = np.arange(k_0, k_f)
    plt.plot(t, y_0, '-k')
    plt.plot(t, y_1, '--r')
    plt.xlabel('Time')
    fig.savefig(loc + 'results_test/flow.png', dpi=300)  # TODO: Naming scheme for files...
    plt.close()


def plot_comparison(y_0: np.ndarray, y_1: np.ndarray, i: List, k_0: int, k_f: int, loc: str):
    fig = plt.figure()
    t = np.arange(k_0, k_f)
    plt.plot(t, y_0[:, i], '-k')
    plt.plot(t, y_1[:, i], '--r')
    plt.xlabel('Time')
    fig.savefig(loc + 'results_test/closed_loop.png', dpi=300)  # TODO: Naming scheme for files...
    plt.close()


def plot_lp(y_0: np.ndarray, i: List, k_0: int, k_f: int, n_update: int):
    fig = plt.figure()
    t = np.arange(k_0, k_f)
    plt.plot(t, y_0[:, i], '-k')
    plt.xlabel('Time')
    plt.ylim((0, 2100))
    fig.savefig('/Users/justshim/src/semester-project/CTM-s/sandbox/figures/results_test/lp' + str(n_update) + '.png', dpi=300)
    plt.close()
