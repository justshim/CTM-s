import os
import numpy as np
from matplotlib import pyplot as plt

from src.results import TrafficPerformance

if __name__ == '__main__':
    """
    Plot results from main.py
    """

    cwd = os.getcwd()
    path = os.path.split(cwd)
    date = "05_12_2023_09_58_57"
    results_path = f"{path[0]}/results/{date}/"
    fig_path = path[0] + '/sandbox/figures/results/'

    k_0 = 2400
    k_f = 4000

    n = 1
    i = 7

    results_0 = np.load(f"{results_path}f_{0}_it_{n}_cell_{i}.npz")
    results_1 = np.load(f"{results_path}f_{1}_it_{n}_cell_{i}.npz")

    # results_0 = np.load(f"{results_path}f_{0}_it_{n}_stat_{0}.npz")
    # results_1 = np.load(f"{results_path}f_{1}_it_{n}_stat_{0}.npz")

    # TODO: Compute demand and supply





    fig = plt.figure()
    t = np.arange(k_0, k_f)

    plt.plot(t, results_0["rho"][k_0:k_f], '-k')
    plt.plot(t, results_1["rho"][k_0:k_f], '--r')
    plt.xlabel('Time')
    fig.savefig(fig_path + 'closed_loop.png', dpi=300)
    plt.close()
