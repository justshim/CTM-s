import itertools
import random

import numpy as np
import pandas as pd
import vaex

if __name__ == '__main__':
    generation = False
    random = True
    # print("GPU test: " + str(tf.test.is_gpu_available()) + "\n")
    path_file_output = "../data/saaaaaaa.csv"
    if (generation == True):
        print("Start: ")
        aa = 16
        bb = aa + 1
        a = range(aa, aa + 1, 1)  # i
        b = range(bb, bb + 4, 1)  # j
        c = range(90, 721, 90)  # delta
        d = range(1, 21, 2)  # beta
        e = range(1, 21, 2)  # priority
        f = range(100, 1001, 200)  # L
        g = range(80, 121, 10)  # v
        h = range(10, 51, 10)  # w
        i = range(1500, 2501, 200)  # qmax
        l = range(50, 101, 10)  # rhomax

        paramlist = list(itertools.product(a, b, c, d, e, f, g, h, i, l))  # all possible combinations of tuples
        print("Total cases to be evaluated: " + str(len(paramlist)))

        arr = np.array(paramlist)
        with open(path_file_output, 'a') as f:
            np.savetxt(f, arr, delimiter=";", fmt='%g')
    if generation is False and random is False:
        # file_path='C:/Users/dspal/Desktop/cazzillo/super_mega_enorme_dataset.csv'
        # dv = vaex.from_csv(file_path, convert=True, chunk_size=10_000_000, sep=';', names = ['i', 'j', 'delta', 'beta', 'priority', 'L_i', 'v_i', 'w_i', 'qmax_i','rhomax_i'])
        dv = vaex.open('C:/Users/User/Desktop/cazzillo/super_mega_enorme_dataset.csv.hdf5')
        xx = np.random.randint(-1, 2, vaex.dataframe.DataFrame.count(dv))
        dv['xx'] = xx
        dv['L_j'] = dv['L_i'] + (dv['L_i'] * 0.25 * dv['xx'])
        dv['v_j'] = dv['v_i'] + (dv['v_i'] * 0.05 * dv['xx'])
        dv['w_j'] = dv['w_i'] + (dv['w_i'] * 0.05 * dv['xx'])
        dv['qmax_j'] = dv['qmax_i'] + (dv['qmax_i'] * 0.10 * dv['xx'])
        dv['rhomax_j'] = dv['rhomax_i'] + (dv['rhomax_i'] * 0.20 * dv['xx'])
        print(dv)
        # dv.export_many('C:/Users/dspal/Desktop/cazzillo/output_chunk-{i:02}.csv', chunk_size=1_000_000)
        # dv.export_hdf5('C:/Users/dspal/Desktop/cazzillo/output_data.hdf5')

    if random is True:
        #dv = vaex.open('C:/Users/User/Desktop/output_data.hdf5')
        dv = vaex.open('C:/Users/User/Desktop/cazzillo/output_data.hdf5')
        values = np.random.randint(0, vaex.dataframe.DataFrame.count(dv), 330000)
        dv2 = dv.take(values)
        print(dv2)
        dv2.export_hdf5('C:/Users/User/Desktop/cazzillo/output_data_random.hdf5')
