from dataclasses import dataclass
from typing import List

import numpy as np
from matplotlib import pyplot as plt

from model.supervisor import Stretch


@dataclass
class TrafficPerformance:
    """
    Class to define, compute performance metrics
    """

    delta_ttt: float        # Additional travel time due to congestion [veh hr]
    delta_ttd: float        # Reduced travel distance due to congestion [veh km]
    e_time: float           # Minutes in the Day with Queue [min]
    e_wait: float           # Average Queue Wait Time [min]


class TrafficEvaluator:
    """
    TODO: ...
    """

    def evaluate(self, s: Stretch) -> TrafficPerformance:
        """
        TODO: ...
        """

        delta_ttt = self.compute_delta_ttt(s)
        delta_ttd = self.compute_delta_ttd(s)
        e_time = self.compute_e_time(s)
        e_wait = self.compute_e_wait(s)

        return TrafficPerformance(delta_ttt, delta_ttd, e_time, e_wait)

    @classmethod
    def compute_delta_ttt(cls, s: Stretch):
        """
        TODO: ...
        """

        delta_ttt = np.zeros_like(s.cells[0].v, dtype=float)

        for i in range(s.n_cells - 1):
            l = s.cells[i].l
            v = np.array(s.cells[i].v, dtype=float)
            v_free = s.cells[i].v_free

            delta_ttt_i = 60 * (np.divide(l, v) - (l / v_free))
            delta_ttt_i[delta_ttt_i < 0] = 0
            delta_ttt += delta_ttt_i

        return np.sum(delta_ttt)

    @classmethod
    def compute_delta_ttd(cls, s: Stretch):
        """
        TODO: ...
        """

        return 0

    @classmethod
    def compute_e_time(cls, s: Stretch):
        """
        TODO: ...
        """

        return 0

    @classmethod
    def compute_e_wait(cls, s: Stretch):
        """
        TODO: ...
        """

        return 0


class TrafficResults:
    """
    TODO: Class to keep track of ...
    """

    n_iter: int
    n_fact: int
    results: List

    def __init__(self, n_iter: int, n_fact: int):
        self.n_iter = n_iter
        self.n_fact = n_fact
        self.results = [[] for _ in range(self.n_fact)]  # TODO: Kinda messy

    def log(self, i: int, p: TrafficPerformance, print_perf: bool):
        """
        TODO: ...
        """

        self.results[i] += [p]
        if print_perf:
            print(p)

    def plot(self):
        """
        TODO: ...
        """
        fig = plt.figure()
        t = np.arange(0, self.n_iter)
        fig.savefig()
        plt.close()


def plot_comparison_test(y_0: np.ndarray, y_1: np.ndarray, i: List, k_0: int, k_f: int, loc: str):
    fig = plt.figure()
    t = np.arange(k_0, k_f)
    plt.plot(t, y_0[:, i], '-k')
    plt.plot(t, y_1[:, i], '--r')
    plt.xlabel('Time')
    fig.savefig(loc + 'comp/cell.png', dpi=300)  # TODO: Naming scheme for files...
    plt.close()


def plot_comparison(y_0: np.ndarray, y_1: np.ndarray, i: List, k_0: int, k_f: int, loc: str):
    fig = plt.figure()
    t = np.arange(k_0, k_f)
    plt.plot(t, y_0[:, i], '-k')
    plt.plot(t, y_1[:, i], '--r')
    plt.xlabel('Time')
    fig.savefig(loc + 'comp/closed_loop.png', dpi=300)  # TODO: Naming scheme for files...
    plt.close()


def plot_lp(y_0: np.ndarray, i: List, k_0: int, k_f: int, n_update: int):
    fig = plt.figure()
    t = np.arange(k_0, k_f)
    plt.plot(t, y_0[:, i], '-k')
    plt.xlabel('Time')
    fig.savefig('/Users/justshim/src/semester-project/CTM-s/sandbox/figures/comp/lp' + str(n_update) + '.png', dpi=300)
    plt.close()
