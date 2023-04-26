import os
from model.factory import Factory
from model.parameters import CTMsParameters
from src.data import TrafficFlowGenerator
from optimizer import TrafficOptimizer
from control import ControlParameters, TrafficHistory
from results import TrafficEvaluator, TrafficResults, plot_comparison


if __name__ == '__main__':
    ####################
    # Input Parameters #
    ####################

    dt = 10                         # Sampling rate [sec]: Default = 10 seconds
    window_length = 90              # Length of one LP window [# of time increments]: 90 = 15 minutes
    freq = 30                       # Frequency of Solving the LP [# of time increments]: 30 = 5 minutes
    max_updates = 20                # Number of control updates: 15
    n_fact = 2                      # Number of factories / algorithms
    n_iter = 2                      # Number of iterations

    day_start = 0                   # Beginning of day [Index in time increments]
    congestion_start = 2400         # Beginning of congestion/control period [Index in time increments]
    congestion_end = 3800           # End of congestion/control period [Index in time increments]
    day_end = 8640                  # Beginning of day [Index in time increments]

    ###################
    # Read Parameters #
    ###################

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

    ########################################
    # Build Factory, LP Solver, Controller #
    ########################################

    fac_0 = Factory(parameters)  # Nominal highway w/o control
    fac_1 = Factory(parameters)  # Highway w/ receding-horizon control
    # fac_2 = Factory(parameters)  # TODO: Highway w/ receding-horizon ILC control

    gen = TrafficFlowGenerator(phi_loc)
    control_parameters = ControlParameters()
    opt = TrafficOptimizer(parameters, control_parameters)
    hist = TrafficHistory(n_iter, window_length, max_updates)
    perf = TrafficEvaluator()
    results = TrafficResults(n_iter, n_fact)

    ##############################
    # Simulate Factory, Solve LP #
    ##############################

    for n in range(n_iter):
        # Different traffic flow (into first cell) for each iteration/day
        phi_day = gen.get_flow(t_0=day_start, t_f=day_end, perturb=False)

        # Initialize Factory Simulation
        fac_0.stretches[0].simulate_init(phi_day)
        fac_1.stretches[0].simulate_init(phi_day)

        # Simulate from day start to beginning of congestion w/o control
        for k in range(day_start, congestion_start):
            fac_0.stretches[0].update(k)
            fac_1.stretches[0].update(k)

        # During congestion period, run control
        for k in range(congestion_start, congestion_end):
            n_update, no_control_update = divmod(k - congestion_start, freq)

            if not no_control_update and n_update < max_updates:
                window_start = k
                window_end = k + window_length

                # No control in first iteration
                if n != 0:
                    # Get initial state
                    x_0 = fac_1.stretches[0].get_state(window_start)
                    s_s_0 = fac_1.stretches[0].get_station_inflow(window_start)

                    # Generate estimate of flow from data
                    phi_0 = hist.estimate_flow()

                    # Solve LP
                    opt.solve_init(x_0, phi_0, s_s_0, window_start, window_length)
                    opt.solve(print_sol=False, plot_sol=True, n_update=n_update)

                    # Update service stations with ramp-meter control term from LP solution
                    fac_1.stretches[0].update_control(window_start, opt.u_rs_c)

                # Log observed vehicle flow into first cell
                hist.log_flow(phi_day[window_start:window_end])

            fac_0.stretches[0].update(k)
            fac_1.stretches[0].update(k)

        # Simulate from end of congestion to end of day w/o control
        for k in range(congestion_end, day_end):
            fac_0.stretches[0].update(k)
            fac_1.stretches[0].update(k)

        # TODO: Track results
        results.log(0, perf.evaluate(fac_0.stretches[0]), print_perf=True)
        results.log(1, perf.evaluate(fac_1.stretches[0]), print_perf=True)

        # Reset factories at end of each iteration
        fac_0.stretches[0].reset()
        fac_1.stretches[0].reset()

    # TODO: Plot final results
    # results.plot()

    # y = fac.stretches[0].y_rho[congestion_start:congestion_end, :]
    # y = fac.stretches[0].u_rs[congestion_start:congestion_end, :]
    # y = fac.stretches[0].y_e[congestion_start:congestion_end, :]

    t_0 = 2400
    t_f = 3200

    y_0 = fac_0.stretches[0].u_rs[t_0:t_f, :]
    y_1 = fac_1.stretches[0].u_rs[t_0:t_f, :]

    # Plotting
    plot_comparison(y_0, y_1, [0], t_0, t_f, fig_path)

