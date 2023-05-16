import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

from src.results import TrafficPerformance

if __name__ == '__main__':
    """
    Plot results from main.py
    """

    cwd = os.getcwd()
    path = os.path.split(cwd)
    # date = "05_16_2023_10_41_30"  # 2% reduction
    # date = "05_16_2023_12_31_00"  # 12% reduction
    date = "05_16_2023_15_12_14"  # 23% reduction

    results_path = f"{path[0]}/results/{date}/"
    fig_path = path[0] + '/sandbox/figures/results/'

    k_0 = 2400
    k_f = 4000
    plot_times = True

    n = 1
    i = 10

    results_0 = np.load(f"{results_path}f_{0}_it_{n}_cells.npz")
    results_0_s = np.load(f"{results_path}f_{0}_it_{n}_stations.npz")
    results_1 = np.load(f"{results_path}f_{1}_it_{n}_cells.npz")
    results_1_s = np.load(f"{results_path}f_{1}_it_{n}_stations.npz")

    perf_0 = TrafficPerformance.from_array(np.load(f"{results_path}f_{0}_it_{n}_perf.npy"))
    perf_1 = TrafficPerformance.from_array(np.load(f"{results_path}f_{1}_it_{n}_perf.npy"))

    print(perf_0)
    print(perf_1)

    reduction = (perf_0.cong_sev - perf_1.cong_sev) / perf_0.cong_sev
    print(reduction)

    fig = plt.figure()

    if plot_times:
        times = np.arange(np.datetime64('1999-09-24'), np.datetime64('1999-09-25'), np.timedelta64(10, 's'))
        t = times[k_0:k_f]
    else:
        t = np.arange(k_0, k_f)

    # plt.plot(t, results_0_s["r_s"][k_0:k_f, 0], '-k', label="Nominal Station Outflow")
    # plt.plot(t, results_1_s["r_s"][k_0:k_f, 0], '--r', label="Controlled Station Outflow")

    # plt.plot(t, results_0["dem"][k_0:k_f, i], '-k', linewidth=3)
    # plt.plot(t, results_0["sup"][k_0:k_f, i], '-b', linewidth=3)
    # plt.plot(t, results_0["phi"][k_0:k_f, i], '--r')
    # plt.plot(t, results_0["phi"][k_0:k_f, i] + results_0_s["r_s"][k_0:k_f, 0], '--r')

    # plt.plot(t, results_1["dem"][k_0:k_f, i], '-k', linewidth=3, label="Demand")
    # plt.plot(t, results_1["sup"][k_0:k_f, i], '-b', linewidth=3, label="Supply")
    # plt.plot(t, results_1["phi"][k_0:k_f, i], '--r', label="Flow")
    # plt.plot(t, results_1["phi"][k_0:k_f, i] + results_1_s["r_s"][k_0:k_f, 0], '--r')

    plt.plot(t, results_0["rho"][k_0:k_f, 6], '-k', label='Nominal Cell 6 Density')
    plt.plot(t, results_1["rho"][k_0:k_f, 6], '--r', label='Controlled Cell 6 Density')
    plt.plot(t, results_0["rho"][k_0:k_f, 7], '-g', label='Nominal Cell 7 Density')
    plt.plot(t, results_1["rho"][k_0:k_f, 7], '--c', label='Controlled Cell 7 Density')
    plt.plot(t, results_0["rho"][k_0:k_f, 8], '-b', label='Nominal Cell 8 Density')
    plt.plot(t, results_1["rho"][k_0:k_f, 8], '--m', label='Controlled Cell 8 Density')
    plt.plot(t, results_0["rho"][k_0:k_f, 9], '-', label='Nominal Cell 9 Density', color='orange')
    plt.plot(t, results_1["rho"][k_0:k_f, 9], '--y', label='Controlled Cell 9 Density')


    # plt.plot(t, results_0_s["r_s"][k_0:k_f, 0], '-k')
    # plt.plot(t, results_1_s["r_s"][k_0:k_f, 0], '--r')

    if plot_times:
        plt.xlabel('Time of Day')
        fig.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    else:
        plt.xlabel('Time Step')

    plt.ylabel('Vehicle Density [veh/km]')
    plt.legend(loc="upper left", fontsize=7)
    fig.savefig(fig_path + 'closed_loop.png', dpi=300)
    plt.close()
