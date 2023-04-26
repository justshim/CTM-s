import os

from src.data import *
from model.factory import Factory
from model.parameters import CTMsParameters
from optimizer import TrafficOptimizer
from results import plot_comparison_test
from control import ControlParameters


if __name__ == '__main__':
    ####################
    # Input Parameters #
    ####################

    dt = 10                     # Sampling rate [sec]
    window_length = 90         # Length of one interation [hrs]
    n_iter = 1                  # Number of ILC iterations
    eta = 1                     # Optimization Parameter: Tradeoff Between TTS, TTD
    start = 0
    congestion_start = 2630
    congestion_end = congestion_start + window_length
    end = 4000

    ############################
    # Read Parameters and Data #
    ############################

    cwd = os.getcwd()
    path = os.path.split(cwd)
    hi_loc = path[0] + '/data/CTM_param_out_A2_15cells.csv'
    onr_loc = path[0] + '/data/onramps.csv'
    offr_loc = path[0] + '/data/offramps.csv'
    st_loc = path[0] + '/data/stations_one.csv'
    phi_loc = path[0] + '/data/phi_1_24h_realsmooth.csv'
    phi_onr_loc = path[0] + '/data/onramps_signal.csv'
    fig_path = path[0] + '/sandbox/figures/'

    parameters = CTMsParameters(hi_loc, onr_loc, offr_loc, st_loc, phi_onr_loc)

    # parameters.update_dt(dt=dt/3600)  # TODO: Default value

    phi_0 = read_phi(phi_loc, t_0=start, t_f=end, ds=dt / 10)

    ########################################
    # Build Factory, LP Solver, Controller #
    ########################################

    fac = Factory(parameters)

    gen = TrafficFlowGenerator(phi_loc)
    control_parameters = ControlParameters()
    opt = TrafficOptimizer(parameters, control_parameters)
    # hist = TrafficHistory(n_iter, window_length, max_updates)
    # perf = TrafficEvaluator()
    # results = TrafficResults(n_iter, n_fact)

    ##############################
    # Simulate Factory, Solve LP #
    ##############################

    day_length = end - start
    plot_ids = list(range(0, 14))

    for n in range(n_iter):
        # CTM-s Dynamics
        fac.stretches[0].simulate_init(phi_0)

        for k in range(start, congestion_start):
            fac.stretches[0].update(k)

        # During congestion period, run control
        x_0 = fac.stretches[0].get_state(k=congestion_start)
        s_s_0 = fac.stretches[0].get_station_inflow(k=congestion_start)
        phi_0 = read_phi(phi_loc, t_0=congestion_start, t_f=congestion_end, ds=dt/10)  # In reality, should be some average

        # Linear Program
        k_0 = congestion_start
        opt.solve_init(x_0, phi_0, s_s_0, k_0, window_length)
        opt.solve(print_sol=True)
        # y_lp = opt.y_rho[:-1, :]
        # y_lp = opt.u_phi
        # y_lp = opt.u_rs_c
        y_lp = opt.y_e[:-1, :]

        for k in range(congestion_start, congestion_end):
            fac.stretches[0].update(k)

        # y = fac.stretches[0].y_rho[congestion_start:congestion_end, :]
        # y = fac.stretches[0].u_phi[congestion_start:congestion_end, :]
        # y = fac.stretches[0].u_rs[congestion_start:congestion_end, :]
        y = fac.stretches[0].y_e[congestion_start:congestion_end, :]

        # Plotting
        plot_comparison_test(y, y_lp, [0], congestion_start, congestion_end, fig_path)
