import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from model import factory as f


if __name__ == '__main__':
    """ Read in parameters """
    # cwd = os.getcwd()
    # datafile = '/data/CTM_param_out_A2_15cells.csv'
    # path = cwd + datafile

    loc = '/Users/justshim/src/semester-project/CTM-s/data/CTM_param_out_A2_15cells.csv'  # TODO: Update to automatic path generation
    parameters = pd.read_csv(loc, sep=';').to_numpy()

    ID = parameters[:, 0]  # Cell ID
    L = parameters[:, 1]  # Cell Length [km]
    v_free = parameters[:, 2]  # Free-flow Speed [km/hr]
    w_cong = parameters[:, 3]  # Congestion Wave Speed [km/hr]
    q_max = parameters[:, 4]  # Maximum Cell Capacity [veh/hr]
    rho_max = parameters[:, 5]  # Maximum Jam Density [veh/km]
    dt = parameters[:, 6]  # dt [10 sec (1 day/3600 sec) = 1/360 day]
    p_ms = 1.00  # TODO: Update this value! 0.93 w/ service station

    phi_loc = '/Users/justshim/src/semester-project/CTM-s/data/phi_1_24h_realsmooth.csv'  # TODO: Update to automatic path generation
    phi_real = []
    with open(phi_loc, encoding='utf8') as fa:
        for line in fa:
            phi_real.append(float(line.strip()))

    """Create factory object"""
    fac = f.Factory()

    # Create highway stretch
    fac.createStretch(time_length=dt[0],
                      last_phi=phi_real,
                      first_d_big=phi_real)  # TODO: Refactor: time_length to dt

    # Add cells to highway stretch
    for i in range(ID.size):
        fac.addCellToStretch(id_stretch=0,
                             length=L[i],
                             v_free=v_free[i],
                             w=w_cong[i],
                             q_max=q_max[i],
                             rho_max=rho_max[i],
                             p=p_ms)

    # TODO: Add on-ramps, off-ramps to the stretch
    # fac.addOnRampToStretch(id_stretch=0,
    #                        d_r=0,
    #                        r_r_max=0,
    #                        p_r=0)
    #
    # fac.addOffRampToStretch(id_stretch=0,
    #                         i=0,
    #                         beta_r=0)

    # Add service station to stretch
    fac.addStationToStretch(id_stretch=0,
                            r_s_max=1500,
                            i=5,
                            j=7,
                            delta=80,
                            beta_s=0.10,
                            p=0.07)

    # TODO: Better way of doing this?
    fac.stretches[0].cells[7].p_ms = 0.93

    """Generate fundamental diagram density vs flow"""
    # Plot mainline demand from preceding cell
    demand_ms = []
    # Plot mainline supply from current cell
    supply_ms = []

    # Plot on-ramp demand
    demand_ss = []

    # Plot mainline supply - on-ramp demand
    supply_ms_alloc = []

    # Plot mainline supply * priority of mainline
    supply_ms_p = []

    # Plot mainline supply - priority of service station
    supply_ss_p = []

    # TODO: How can I accurately recreate the fundamental diagram?
    #  My first idea was to manually change the density, but in doing,
    #  I break conservation equations...

    """
    # Need to get to a steady state in the cell behavior...
    # Simulate up until around 9 AM, around the first congestion peak
    k_ = 100
    for k in range(k_):  # 9 AM = 3240
        fac.stretches[0].update(k)

    fac_9 = copy.deepcopy(fac)

    # Iterate over the cells
    # Note: Only cell 13 (14?) will have a service station dynamic effect
    # while i in range(1, ID.size):
    i = 13
    rhos = np.linspace(0, fac.stretches[0].cells[i].rho_max, 1000).tolist()
    phis = []

    for rho in rhos:
        fac.stretches[0].cells[i].rho[-1] = rho  # Manually set the density
        fac.stretches[0].update(kappa=k_)  # Update the highway stretch
        print(fac.stretches[0].cells[i].rho[-1])
        phis.append(fac.stretches[0].cells[i].phi)
        demand_ms.append(fac.stretches[0].computeDPrec(i+1))
        supply_ms.append(fac.stretches[0].cells[i].s_big)
        demand_ss.append(fac.stretches[0].stations[0].d_s_big)
        supply_ms_alloc.append(supply_ms[-1] - demand_ss[-1])
        supply_ms_p.append(supply_ms[-1] * fac.stretches[0].cells[i].p_ms)
        supply_ss_p.append(supply_ms[-1] * (1-fac.stretches[0].cells[i].p_ms))

        fac = fac_9  # TODO: Reset mechanism via deep copy? Does this work?
    """

    i = 7
    rhos = []
    phis = []

    for k in range(8640):  # 24 hrs = 8640 "deciseconds"
        fac.stretches[0].update(k)
        rhos.append(fac.stretches[0].cells[i].rho[-1])
        phis.append(fac.stretches[0].cells[i].phi)
        demand_ms.append(fac.stretches[0].computeDPrec(i))
        supply_ms.append(fac.stretches[0].cells[i].s_big)
        demand_ss.append(fac.stretches[0].stations[0].d_s_big)
        supply_ms_alloc.append(supply_ms[-1] - demand_ss[-1])
        supply_ms_p.append(supply_ms[-1] * fac.stretches[0].cells[i].p_ms)
        supply_ss_p.append(supply_ms[-1] * (1 - fac.stretches[0].cells[i].p_ms))

    # Sorting
    order = np.argsort(rhos)
    rhos = np.array(rhos)[order]
    phis = np.array(phis)[order]
    demand_ms = np.array(demand_ms)[order]
    supply_ms = np.array(supply_ms)[order]
    demand_ss = np.array(demand_ss)[order]
    supply_ms_alloc = np.array(supply_ms_alloc)[order]
    supply_ms_p = np.array(supply_ms_p)[order]
    supply_ss_p = np.array(supply_ss_p)[order]


    # Plotting
    fig = plt.figure()
    # plt.plot(rhos, phis)
    plt.plot(rhos, demand_ms, linewidth=2, linestyle='solid', color='orange')
    plt.plot(rhos, supply_ms, linewidth=2, linestyle='solid', color='green')
    plt.plot(rhos, demand_ss, linewidth=2, linestyle='solid', color='blue')

    plt.plot(rhos, supply_ms_alloc, linewidth=0.2, linestyle='solid', color='lime')
    plt.plot(rhos, supply_ms_p, linewidth=0.2, linestyle='solid', color='cyan')
    plt.plot(rhos, supply_ss_p, linewidth=0.2, linestyle='solid', color='lime')

    plt.plot(rhos, phis, linewidth=0.1, linestyle='solid', color='red')


    fig.savefig('/Users/justshim/src/semester-project/CTM-s/figures/fundamental_diagram.png', dpi=300)
