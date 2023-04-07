import os

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

from read_data import read_parameters, read_phi
from model.factory import Factory


if __name__ == '__main__':
    ####################
    # Input Parameters #
    ####################

    start = 0
    end = 8640

    ############################
    # Read Parameters and Data #
    ############################

    cwd = os.getcwd()
    path = os.path.split(cwd)
    hi_loc = path[0] + '/data/CTM_param_out_A2_15cells.csv'
    onr_loc = path[0] + '/data/onramps.csv'
    offr_loc = path[0] + '/data/offramps.csv'
    stations_loc = path[0] + '/data/stations_one.csv'
    phi_loc = path[0] + '/data/phi_1_24h_realsmooth.csv'
    phi_onr_loc = path[0] + '/data/onramps_signal.csv'
    fig_path = path[0] + '/sandbox/figures/ctms_dynamics/realsmooth/stations_one/'

    parameters = read_parameters(hi_loc,
                                 onr_loc,
                                 offr_loc,
                                 stations_loc,
                                 phi_onr_loc)

    phi_0 = read_phi(phi_loc, t_0=start, t_f=end)

    ########################################
    # Build Factory, LP Solver, Controller #
    ########################################

    fac = Factory(parameters)

    ####################
    # Simulate Factory #
    ####################

    fac.stretches[0].simulate_init(phi_0)

    for k in range(start, end):
        fac.stretches[0].update(k)

    #######################
    # Plot CTM-s Dynamics #
    #######################

    # Cell Densities
    fig_1 = plt.figure()
    times = np.arange(np.datetime64('1999-09-24'), np.datetime64('1999-09-25'), np.timedelta64(10, 's'))
    times = times[start:end]
    # times = np.arange(start, end)
    rho = fac.stretches[0].x_rho[:-1, :-1]
    plt.plot(times, rho)
    plt.xlabel('Time')
    plt.ylabel('Vehicle Density [veh/km]')
    fig_1.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig_1.autofmt_xdate()

    fig_1.savefig(fig_path + '/densities.png', dpi=300)

    # Flows
    fig_2 = plt.figure()
    phi = fac.stretches[0].u_phi[:, :-1]
    plt.plot(times, phi)
    plt.xlabel('Time')
    plt.ylabel('Vehicle Flow [veh/hr]')
    fig_2.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig_2.autofmt_xdate()

    fig_2.savefig(fig_path + '/flows.png', dpi=300)

    # Fundamental Diagram
    fig_3 = plt.figure()
    for i in range(0, 14):
        plt.scatter(rho[:, i], phi[:, i], s=1)
        plt.xlabel('Vehicle Density [veh/km]')
        plt.ylabel('Vehicle Flow [veh/hr]')

    fig_3.savefig(fig_path + '/fundamental_diagram.png', dpi=300)

    # Service Station
    fig_4 = plt.figure()
    l = fac.stretches[0].x_l[:-1, :]
    e = fac.stretches[0].x_e[:-1, :]
    plt.plot(times, l)
    plt.plot(times, e)
    plt.xlabel('Time')
    plt.ylabel('Number of Vehicles [veh]')
    fig_4.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig_4.autofmt_xdate()

    fig_4.savefig(fig_path + '/service_station.png', dpi=300)
