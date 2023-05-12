import os

from src.data import *
from model.factory import Factory
from model.parameters import CTMsParameters
from optimizer import TrafficOptimizer
from results import plot_comparison_test, plot_flow_comparison
from control import ControlParameters


if __name__ == '__main__':
    ####################
    # Input Parameters #
    ####################

    dt = 10                     # Sampling rate [sec]
    window_length = 270         # Length of one interation [hrs]
    n_iter = 1                  # Number of ILC iterations
    eta = 1                     # Optimization Parameter: Tradeoff Between TTS, TTD

    day_start = 0
    congestion_start = 2400
    congestion_end = 4000
    day_end = 4000

    cwd = os.getcwd()
    path = os.path.split(cwd)
    hi_loc = f"{path[0]}/data/CTM_param_out_A2_15cells.csv"
    onr_loc = f"{path[0]}/data/onramps.csv"
    offr_loc = f"{path[0]}/data/offramps.csv"
    st_loc = f"{path[0]}/data/stations_one.csv"
    phi_loc = f"{path[0]}/data/phi_1_24h_realsmooth.csv"
    phi_onr_loc = f"{path[0]}/data/onramps_signal.csv"
    fig_path = f"{path[0]}/sandbox/figures/"

    ########################################
    # Build Factory, LP Solver, Controller #
    ########################################

    parameters = CTMsParameters(hi_loc, onr_loc, offr_loc, st_loc, phi_onr_loc)
    fac = Factory(parameters)
    gen = TrafficFlowGenerator(phi_loc, congestion_start, congestion_end)
    control_parameters = ControlParameters()
    opt = TrafficOptimizer(parameters, control_parameters)

    ##############################
    # Simulate Factory, Solve LP #
    ##############################

    day_length = day_end - day_start
    plot_ids = list(range(0, 14))

    phi_day = 1.2 * gen.get_flow(t_0=day_start, t_f=day_end, perturb=False)

    window_start = 2900
    window_end = window_start + window_length

    for n in range(n_iter):
        # CTM-s Dynamics
        fac.stretches[0].simulate_init(phi_day)

        for k in range(day_start, window_start):
            fac.stretches[0].update(k)

        # During congestion period, run control
        x_0 = fac.stretches[0].get_state(k=window_start)
        s_s_0 = fac.stretches[0].get_station_inflow(k=window_start)
        phi_0 = phi_day[window_start:window_end]

        # phi_0 = gen.get_flow(t_0=window_start, t_f=window_end, perturb=True)
        # plot_flow_comparison(phi_day[window_start:window_end], phi_0, window_start, window_end, fig_path)

        # Linear Program
        k_0 = window_start
        opt.solve_init(x_0, phi_0, s_s_0, k_0, window_length)
        opt.solve(print_sol=False)
        # y_lp = opt.y_rho[:-1, :]
        # y_lp = opt.u_phi
        y_lp = opt.u_rs_c
        # y_lp = opt.y_e[:-1, :]

        for k in range(window_start, window_end):
            fac.stretches[0].update(k)

        # y = fac.stretches[0].y_rho[window_start:window_end, :]
        # y = fac.stretches[0].u_phi[window_start:window_end, :]
        y = fac.stretches[0].u_rs[window_start:window_end, :]
        # y = fac.stretches[0].y_e[window_start:window_end, :]

        # Plotting
        plot_comparison_test(y, y_lp, [0], window_start, window_end, fig_path)
