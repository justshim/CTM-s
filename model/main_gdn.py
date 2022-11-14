from model import factory as f
import xlrd


def ga(stretch, phi, duration, rsmax, station):

    phi_zero = phi
    last_phi = phi

    # create a factory instance that manages the creation of objects
    fac0 = f.Factory()
    fac0.createStretch(1/360, last_phi, phi_zero)
    for row in range(stretch.shape[0]):
        fac0.addCellToStretch(0, stretch[row][1], stretch[row][2], stretch[row][3], stretch[row][4], stretch[row][5], 1)

    max_delta0 = first_iteration(fac0, duration)
    #print("max delta0")
    #print(max_delta0)

    fac1 = f.Factory()
    fac1.createStretch(1/360, last_phi, phi_zero)
    for row in range(stretch.shape[0]):
        fac1.addCellToStretch(0, stretch[row][1], stretch[row][2], stretch[row][3], stretch[row][4], stretch[row][5], 1)

    fac1.stretches[0].cells[station[0]+2].p_ms = fac1.stretches[0].cells[station[0]+2].p_ms - 0.05
    fac1.addStationToStretch(0, rsmax, station[0], station[0]+2, station[1], station[2], 0.05)

    result = iteration(fac1, duration)
    #print(result)
    integral = result[0]
    max_delta = result[1]
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
