import os
from datetime import datetime
from model.factory import Factory
from model.parameters import CTMsParameters
from src.data import TrafficFlowGenerator
from optimizer import TrafficOptimizer
from control import ControlParameters, TrafficHistory
from results import TrafficEvaluator, TrafficResults


if __name__ == '__main__':
    ####################
    # Input Parameters #
    ####################

    dt = 10                         # Sampling rate [sec]: Default = 10 seconds
    window_length = 270             # Length of one LP window [# of time  increments]: 270 = 45 minutes
    freq = 60                       # Frequency of Solving the LP [# of time increments]: 60 = 10 minutes
    n_fact = 2                      # Number of factories / algorithms
    n_iter = 2                      # Number of iterations

    day_start = 0                   # Beginning of day [Index in time increments]
    congestion_start = 2400         # Beginning of congestion/control period [Index in time increments]
    congestion_end = 4000           # End of congestion/control period [Index in time increments]
    day_end = 8640                  # Beginning of day [Index in time increments]

    #################
    # Preliminaries #
    #################

    cwd = os.getcwd()
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")

    path = os.path.split(cwd)
    hi_loc = f"{path[0]}/data/CTM_param_out_A2_15cells.csv"
    onr_loc = f"{path[0]}/data/onramps.csv"
    offr_loc = f"{path[0]}/data/offramps.csv"
    st_loc = f"{path[0]}/data/stations_one.csv"
    phi_loc = f"{path[0]}/data/phi_1_24h_realsmooth.csv"
    phi_onr_loc = f"{path[0]}/data/onramps_signal.csv"
    fig_path = f"{path[0]}/sandbox/figures/"

    results_path = f"{path[0]}/sandbox/results/{date_time}/"
    os.mkdir(results_path)

    day_length = day_end - day_start
    max_updates, _ = divmod(congestion_end - congestion_start, freq)

    ######################
    # Initialize Classes #
    ######################

    parameters = CTMsParameters(hi_loc, onr_loc, offr_loc, st_loc, phi_onr_loc)

    fac_0 = Factory(parameters)  # Nominal highway w/o control
    fac_1 = Factory(parameters)  # Highway w/ receding-horizon control (MPC)
    # fac_2 = Factory(parameters)  # TODO: Highway w/ receding-horizon ILC control

    control_parameters = ControlParameters()  # TODO: Highway w/ receding-horizon ILC control
    opt_1 = TrafficOptimizer(parameters, control_parameters)
    # opt_2 = TrafficOptimizer(parameters, control_parameters)  # TODO: Highway w/ receding-horizon ILC control

    gen = TrafficFlowGenerator(phi_loc, congestion_start, congestion_end)
    hist = TrafficHistory(n_iter, window_length, max_updates)
    perf = TrafficEvaluator()
    results = TrafficResults(n_iter, n_fact)

    ##############################
    # Simulate Highway Stretches #
    ##############################

    for n in range(n_iter):
        # Different traffic flow (into first cell) for each iteration/day
        phi_day = 1.2 * gen.get_flow(t_0=day_start, t_f=day_end, perturb=False)

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

                real_phi_0 = phi_day[window_start:window_end]

                # No control in first iteration
                if n != 0:
                    # Get initial state
                    x_0 = fac_1.stretches[0].get_state(window_start)
                    s_s_0 = fac_1.stretches[0].get_station_inflow(window_start)

                    # Generate estimate of flow from data
                    phi_0 = hist.get_flow(predict=False)
                    # plot_flow_comparison(real_phi_0, phi_0, window_start, window_end, fig_path)

                    # Solve LP
                    opt_1.solve_init(x_0, phi_0, s_s_0, window_start, window_length)
                    opt_1.solve(print_sol=False, plot_sol=False, n_update=n_update)

                    # Update service stations with ramp-meter control term from LP solution
                    fac_1.stretches[0].update_control(window_start, opt_1.u_rs_c)

                # Log observed vehicle flow into first cell
                hist.log_flow(real_phi_0)

            fac_0.stretches[0].update(k)
            fac_1.stretches[0].update(k)

        # Simulate from end of congestion to end of day w/o control
        for k in range(congestion_end, day_end):
            fac_0.stretches[0].update(k)
            fac_1.stretches[0].update(k)

        # Calculate performance results
        perf_0 = perf.evaluate(fac_0.stretches[0], print_perf=True)
        perf_1 = perf.evaluate(fac_1.stretches[0], print_perf=True)

        # Save simulation data, performance results
        results.save(0, n, fac_0.stretches[0], perf_0, results_path)
        results.save(1, n, fac_1.stretches[0], perf_1, results_path)

        # Reset factories at end of each iteration
        fac_0.stretches[0].reset()
        fac_1.stretches[0].reset()
