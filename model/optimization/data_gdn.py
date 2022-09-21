import numpy as np

N_STRETCHES = 100000
inputs = []
with open('output.csv', "ab") as f:
    for index in range(N_STRETCHES):
        N_CELLS = np.random.randint(7, 16)
        stretch = np.zeros([N_CELLS, 6])
        for row in range(N_CELLS):
            stretch[row, 0] = row + 1
            stretch[row, 1] = np.random.randint(300, 1001)  # L
            stretch[row, 2] = np.random.randint(80, 121)  # v
            stretch[row, 3] = np.random.randint(10, 41)  # w
            stretch[row, 4] = np.random.randint(1500, 2501)  # q_max
            stretch[row, 5] = np.random.randint(60, 101)  # rho_max
        np.savetxt(f, stretch, fmt='%g', delimiter=';', newline='\n')







