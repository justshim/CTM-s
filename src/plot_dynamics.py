import os

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

from src.data import TrafficFlowGenerator
from model.parameters import CTMsParameters
from model.factory import Factory

if __name__ == '__main__':
    """
    Generate nominal CTM-s Dynamics figures
    """

    ####################
    # Input Parameters #
    ####################

    day_start = 0                   # Beginning of day [Index in time increments]
    congestion_start = 2400         # Beginning of congestion/control period [Index in time increments]
    congestion_end = 3800           # End of congestion/control period [Index in time increments]
    day_end = 8640                  # Beginning of day [Index in time increments]

    with_stations = False

    cwd = os.getcwd()
    path = os.path.split(cwd)
    hi_loc = path[0] + '/data/CTM_param_out_A2_15cells.csv'
    onr_loc = path[0] + '/data/onramps.csv'
    offr_loc = path[0] + '/data/offramps.csv'
    st_loc = path[0] + '/data/stations_none.csv'
    phi_loc = path[0] + '/data/phi_1_24h_realsmooth.csv'
    phi_onr_loc = path[0] + '/data/onramps_signal.csv'
    fig_path = path[0] + '/figures/dynamics/realsmooth/stations_none'

    ########################################
    # Build Factory, LP Solver, Controller #
    ########################################

    parameters = CTMsParameters(hi_loc, onr_loc, offr_loc, st_loc, phi_onr_loc)
    fac = Factory(parameters)
    gen = TrafficFlowGenerator(phi_loc, congestion_start, congestion_end)

    ####################
    # Simulate Factory #
    ####################

    phi_day = gen.get_flow(t_0=day_start, t_f=day_end, perturb=False)

    fac.stretches[0].simulate_init(phi_day)

    for k in range(day_start, day_end):
        fac.stretches[0].update(k)

    #######################
    # Plot CTM-s Dynamics #
    #######################

    # times = np.arange(day_start, day_end)
    times = np.arange(np.datetime64('1999-09-24'), np.datetime64('1999-09-25'), np.timedelta64(10, 's'))
    times_morning = times[congestion_start:congestion_end]

    # Plot Raw Traffic Flow Profile
    fig = plt.figure()
    plt.plot(times, phi_day)
    plt.xlabel('Time of Day')
    plt.ylabel('Vehicle Flow [veh/hr]')
    fig.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig.autofmt_xdate()
    fig.savefig(fig_path + '/traffic_flow_profile.png', dpi=300)
    plt.clf()

    # Fundamental Diagram
    rho = fac.stretches[0].y_rho[:-1, :-1]
    phi = fac.stretches[0].u_phi[:, :-1]

    for i in range(0, 14):
        plt.scatter(rho[:, i], phi[:, i], s=1)
        plt.xlabel('Vehicle Density [veh/km]')
        plt.ylabel('Vehicle Flow [veh/hr]')

    fig.savefig(fig_path + '/fundamental_diagram.png', dpi=300)
    plt.clf()

    # Cell Densities (Full Day)
    plt.plot(times, rho)
    plt.xlabel('Time of Day')
    plt.ylabel('Vehicle Density [veh/km]')
    fig.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig.autofmt_xdate()
    fig.savefig(fig_path + '/cell_densities.png', dpi=300)
    plt.clf()

    # Cell Flows (Full Day)
    plt.plot(times, phi)
    plt.xlabel('Time of Day')
    plt.ylabel('Vehicle Flow [veh/hr]')
    fig.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig.autofmt_xdate()
    fig.savefig(fig_path + '/cell_flows.png', dpi=300)
    plt.clf()

    # Cell Densities (Morning Congestion Period)
    rho = fac.stretches[0].y_rho[congestion_start:congestion_end, :-1]
    plt.plot(times_morning, rho)
    plt.xlabel('Time of Day')
    plt.ylabel('Vehicle Density [veh/km]')
    fig.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig.autofmt_xdate()
    fig.savefig(fig_path + '/cell_densities_morning.png', dpi=300)
    plt.clf()

    # Cell Flows (Morning Congestion Period)
    phi = fac.stretches[0].u_phi[congestion_start:congestion_end, :-1]
    plt.plot(times_morning, phi)
    plt.xlabel('Time of Day')
    plt.ylabel('Vehicle Flow [veh/hr]')
    fig.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig.autofmt_xdate()
    fig.savefig(fig_path + '/cell_flows_morning.png', dpi=300)
    plt.clf()

    if with_stations:
        # Service Station (Full Day)
        l = fac.stretches[0].y_l[:-1, :]
        e = fac.stretches[0].y_e[:-1, :]
        plt.plot(times, l)
        plt.plot(times, e)
        plt.xlabel('Time of Day')
        plt.ylabel('Number of Vehicles [veh]')
        fig.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        fig.autofmt_xdate()
        fig.savefig(fig_path + '/station_state.png', dpi=300)
        plt.clf()

        # Service Station Outflow (Full Day)
        r_s = fac.stretches[0].u_rs[:, 0]
        plt.plot(times, r_s)
        plt.xlabel('Time of Day')
        plt.ylabel('Vehicle Flow [veh/hr]')
        fig.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        fig.autofmt_xdate()
        fig.savefig(fig_path + '/station_flow.png', dpi=300)
        plt.clf()

        # Service Station State (Morning Congestion)
        l = fac.stretches[0].y_l[congestion_start:congestion_end, :]
        e = fac.stretches[0].y_e[congestion_start:congestion_end, :]
        plt.plot(times_morning, l)
        plt.plot(times_morning, e)
        plt.xlabel('Time of Day')
        plt.ylabel('Number of Vehicles [veh]')
        fig.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        fig.autofmt_xdate()
        fig.savefig(fig_path + '/station_state_morning.png', dpi=300)
        plt.clf()

        # Service Station Outflow (Morning Congestion)
        r_s = fac.stretches[0].u_rs[congestion_start:congestion_end, 0]
        plt.plot(times_morning, r_s)
        plt.xlabel('Time of Day')
        plt.ylabel('Vehicle Flow [veh/hr]')
        fig.axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        fig.autofmt_xdate()
        fig.savefig(fig_path + '/station_flow_morning.png', dpi=300)
        plt.clf()

    plt.close('all')
