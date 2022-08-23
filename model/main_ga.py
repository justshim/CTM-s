from model import factory as f
import xlrd


def ga(path, duration, station, onramps, offramps):
    #station = [rsmax, ii, jj, delta, beta, p]

    # File read
    wb = xlrd.open_workbook(path)
    sh = wb.sheet_by_name("Cells parameters")
    sh.cell_value(0, 0)

    # read phi first cell from xls file
    wb_phi = xlrd.open_workbook(path)

    sh_phi = wb_phi.sheet_by_name("First Demand Real")
    sh_phi.cell_value(0, 0)
    phi_zero = []

    for i in range(0, sh_phi.nrows):
        phi_zero.append(sh_phi.cell_value(i, 0))

    sh_last_phi = wb.sheet_by_name("Last Demand Real")
    sh_last_phi.cell_value(0, 0)
    last_phi = []
    for i in range(0, sh_last_phi.nrows):
        last_phi.append(sh_last_phi.cell_value(i, 0))

    # create a factory instance that manages the creation of objects
    fac0 = f.Factory()
    fac0.createStretch(sh.cell_value(2, 6), last_phi, phi_zero)
    for i in range(1, sh.nrows):
        fac0.addCellToStretch(0, sh.cell_value(i, 1), sh.cell_value(i, 2), sh.cell_value(i, 3), sh.cell_value(i, 4), sh.cell_value(i, 5), 1)

    if len(onramps) > 0:
        # create the on-ramps via the factory
        # ID stretch, d_r, r_r_max, j, p_r
        for onr in onramps:
            fac0.addOnRampToStretch(0, onr[0], onr[1], onr[2], onr[3])
            fac0.stretches[0].cells[onr[2]].p_ms = fac0.stretches[0].cells[onr[2]].p_ms - onr[3]

    if len(offramps) > 0:
        # create the off-ramps via the factory
        # ID_stretch, i, beta_r
        for offr in offramps:
            fac0.addOffRampToStretch(0, offr[0], offr[1])

    max_delta0 = first_iteration(fac0, duration)

    fac1 = f.Factory()
    fac1.createStretch(sh.cell_value(2, 6), last_phi, phi_zero)
    for i in range(1, sh.nrows):
        fac1.addCellToStretch(0, sh.cell_value(i, 1), sh.cell_value(i, 2), sh.cell_value(i, 3),
                              sh.cell_value(i, 4), sh.cell_value(i, 5), 1)

    fac1.stretches[0].cells[station[1]].p_ms = fac1.stretches[0].cells[station[1]].p_ms - station[5]
    fac1.addStationToStretch(0, station[0], station[1], station[2], station[3], station[4], station[5])

    if len(onramps) > 0:
        # create the on-ramps via the factory
        # ID stretch, d_r, r_r_max, j, p_r
        for onr in onramps:
            fac1.addOnRampToStretch(0, onr[0], onr[1], onr[2], onr[3])
            fac1.stretches[0].cells[onr[2]].p_ms = fac1.stretches[0].cells[onr[2]].p_ms - onr[3]

    if len(offramps) > 0:
        # create the off-ramps via the factory
        # ID_stretch, i, beta_r
        for offr in offramps:
            fac1.addOffRampToStretch(0, offr[0], offr[1])

    result = iteration(fac1, duration)
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
