import os
from read_data import read_parameters, read_phi
from model.factory import Factory
from optimizer import MyTrafficOptimizer


if __name__ == '__main__':
    ####################
    # Input Parameters #
    ####################

    dt = 10                         # Sampling rate [sec]
    window_length = 180             # Length of one LP window [# TODO: hrs]
    freq = 30                       # Frequency of Solving the LP []
    n_iter = 1                      # Number of ILC iterations
    eta = 1                         # Optimization Parameter: Tradeoff Between TTS, TTD

    day_start = 0                   # TODO
    congestion_start = 2160         # TODO
    congestion_end = 4000           # TODO
    day_end = 4000                  # TODO

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

    parameters.update_dt(dt=dt / 3600)

    ########################################
    # Build Factory, LP Solver, Controller #
    ########################################

    fac = Factory(parameters)
    opt = MyTrafficOptimizer(parameters, eta)  # TODO: Eta + Other LP parameters! Make a Control Parameter Class...
    # cont = ILCTrafficController()

    ##############################
    # Simulate Factory, Solve LP #
    ##############################

    ds = int(dt/10)
    phi_day = read_phi(phi_loc, t_0=day_start, t_f=day_end, ds=ds)
    # day_length = day_end - day_start
    # results = MyTrafficResults(day_length)

    for n in range(n_iter):
        # Initialize Factory Simulation
        fac.stretches[0].simulate_init(phi_day)

        # Simulate from day start to beginning of congestion w/o control
        for k in range(day_start, congestion_start):
            fac.stretches[0].update(k)

        # During congestion period, run control
        for k in range(congestion_start, congestion_end):
            (j, no_control_update) = divmod(k - congestion_start, freq)
            if not no_control_update:
                window_start = k
                window_end = k + window_length
                x_0 = fac.stretches[0].get_state(k)
                phi_s_0 = fac.stretches[0].get_station_inflow(k=congestion_start)
                phi_0 = read_phi(phi_loc, t_0=window_start, t_f=window_end, ds=ds)

                opt.solve_init(x_0, phi_0, phi_s_0, window_start, window_length)  # TODO: Check if k_0 is correct..
                opt.solve(print_sol=False)
                # TODO: Keep track of results

                u_rs_c = opt.get_control()  # TODO?
                fac.stretches[0].update_control(u_rs_c)  # TODO?

            fac.stretches[0].update(k)
            # TODO: Keep Track of Results

        # Simulate from end of congestion to end of day w/o control
        for k in range(congestion_end, day_end):
            fac.stretches[0].update(k)
            # TODO: Keep Track of Results

        # TODO: Results, Plotting

        # Reset Factory
        fac.stretches[0].reset()
