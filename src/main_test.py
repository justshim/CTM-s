import os

from read_data import read_parameters, read_phi
from model.factory import Factory
from optimizer import MyTrafficOptimizer
from results import plot_comparison


if __name__ == '__main__':
    ####################
    # Input Parameters #
    ####################

    dt = 10                     # Sampling rate [sec]
    window_length = 100         # Length of one interation [hrs]
    n_iter = 1                  # Number of ILC iterations
    eta = 1                     # Optimization Parameter: Tradeoff Between TTS, TTD
    start = 0
    congestion_start = 2500
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
    stations_loc = path[0] + '/data/stations_one.csv'
    phi_loc = path[0] + '/data/phi_1_24h_realsmooth.csv'
    phi_onr_loc = path[0] + '/data/onramps_signal.csv'
    fig_path = path[0] + '/sandbox/figures/'

    parameters = read_parameters(hi_loc,
                                 onr_loc,
                                 offr_loc,
                                 stations_loc,
                                 phi_onr_loc)

    # parameters.update_dt(dt=dt/3600)  # TODO: Default value

    phi_0 = read_phi(phi_loc, t_0=start, t_f=end, ds=dt / 10)

    ########################################
    # Build Factory, LP Solver, Controller #
    ########################################

    fac = Factory(parameters)
    opt = MyTrafficOptimizer(parameters, eta)

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
        phi_s_0 = fac.stretches[0].get_station_inflow(k=congestion_start)
        phi_0 = read_phi(phi_loc, t_0=congestion_start, t_f=congestion_end, ds=dt / 10)

        # Linear Program
        k_0 = congestion_start
        opt.solve_init(x_0, phi_0, phi_s_0, k_0, window_length)
        opt.solve(print_sol=True)
        y_lp = opt.x_rho[:-1, :]
        # y_lp = opt.u_phi

        for k in range(congestion_start, congestion_end):
            fac.stretches[0].update(k)

        y = fac.stretches[0].x_rho[congestion_start:congestion_end, :]
        # y = fac.stretches[0].u_phi[congestion_start:congestion_end, :]

        # Plotting
        plot_comparison(y, y_lp, plot_ids, congestion_start, congestion_end, fig_path)
