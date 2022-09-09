from model import factory as f
import xlrd
import numpy as np
import random

N_CELLS = 20
DURATION = 8640
def simulation(params, path_phi):
    i = int(params[0])
    j = int(params[1])

    delta = int(params[2])
    beta = int(params[3])/100
    priority_s = int(params[4])/100
    #print(priority_s)

    L_i = params[5]
    v_i = int(params[6])
    w_i = int(params[7])
    qmax_i = int(params[8])
    rhomax_i = int(params[9])

    L_j = params[11]
    v_j = params[12]
    w_j = params[13]
    qmax_j = params[14]
    rhomax_j = params[15]

    # read phi first cell from xls file
    wb_phi = xlrd.open_workbook(path_phi)

    sh_phi = wb_phi.sheet_by_name("First Demand Real")
    sh_phi.cell_value(0, 0)
    phi_zero = []

    for t in range(0, sh_phi.nrows):
        phi_zero.append(sh_phi.cell_value(t, 0))

    sh_last_phi = wb_phi.sheet_by_name("Last Demand Real")
    sh_last_phi.cell_value(0, 0)
    last_phi = []
    for t in range(0, sh_last_phi.nrows):
        last_phi.append(sh_last_phi.cell_value(t, 0))

    # create a factory instance that manages the creation of objects
    fac0 = f.Factory()
    fac0.createStretch(0.00277777777777778, last_phi, phi_zero)

    length = np.random.randint(min(L_i, L_j)*1.2, max(L_i, L_j)*1.6, 20) / 1000
    length[i] = L_i/1000
    length[j] = L_j/1000
    v_free = np.random.randint((v_i - 0.15 * v_i), (v_i + 0.15 * v_i), 20)
    w = np.random.randint((w_i - 0.15 * w_i), (w_i + 0.15 * w_i), 20)
    q_max = np.random.randint((qmax_i - 0.15 * qmax_i), (qmax_i + 0.15 * qmax_i), 20)
    rho_max = np.random.randint((rhomax_i - 0.15 * rhomax_i), (rhomax_i + 0.15 * rhomax_i), 20)
    #print(length)
    #length = np.full((20, 1), 1)

    for t in range(N_CELLS):
        #print(length[t], v_free[t], w[t], q_max[t], rho_max[t])
        fac0.addCellToStretch(0, length[t], v_free[t], w[t], q_max[t], rho_max[t], 1)

    max_delta0 = first_iteration(fac0, DURATION)

    fac1 = f.Factory()
    fac1.createStretch(0.00277777777777778, last_phi, phi_zero)

    for t in range(N_CELLS):
        fac1.addCellToStretch(0, length[t], v_free[t], w[t], q_max[t], rho_max[t], 1)

    fac1.stretches[0].cells[j].p_ms = fac1.stretches[0].cells[j].p_ms - priority_s
    fac1.addStationToStretch(0, 1000, i, j, delta, beta, priority_s)

    result = iteration(fac1, DURATION)
    integral = result[0]
    max_delta = result[1]

    if(max_delta0 == 0):
        print("No congestion")
        pi = 1
    else:
        pi = (max_delta0 - max_delta) / max_delta0


    return integral, pi


def iteration(fac1, duration):
    fac1.stretches[0].computeTTT()
    d = []

    k = 0
    while k < duration:  # k=24h=8640 , k=1h=360, k=3h=1080
        fac1.stretches[0].update(k)
        k = k + 1
        d.append(fac1.stretches[0].delta_big[k - 1])

    integral = 0
    for delta in d:
        integral = integral + delta

    delta_max = max(d)
    return integral, delta_max


def first_iteration(fac0, duration):
    fac0.stretches[0].computeTTT()
    d0 = []
    k = 0
    while k < duration:  # execution of the first simulation without any station
        fac0.stretches[0].update(k)
        k = k + 1
        d0.append(fac0.stretches[0].delta_big[k - 1])

    delta_max0 = max(d0)
    return delta_max0
