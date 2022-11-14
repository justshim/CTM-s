import csv
import itertools
import random
from math import ceil
import time
import numpy as np
import pandas as pd
import vaex
import threading
import multiprocessing
import os

from p_tqdm import p_map

from model.optimization import main_mega as mm


if __name__ == '__main__':
    dv = vaex.open('C:/Users/User/Desktop/cazzillo/output_data_random.hdf5')
    #path_phi = 'C:/Users/User/Documents/MATLAB/CTMs-identification/fnc/extracted_data/CTM_param_out_nice.xls'
    t = time.time()
    # row = dv.take([1])
    # a = row.values.flatten()
    # print(row)
    integral = np.zeros(vaex.dataframe.DataFrame.count(dv))
    pi = np.zeros(vaex.dataframe.DataFrame.count(dv))
    #print(results)
    #for idx in range(vaex.dataframe.DataFrame.count(dv)):
   #     row = dv.take([idx])
     #   print(idx)
       # integral[idx], pi[idx] = mm.simulation(row.values.flatten())
    with open("C:/Users/User/Desktop/cazzillo/output_sim1.csv", 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f, delimiter=';')
        header = ["i", "j", "delta", "beta", "priority", "L_i", "v_i", "w_i", "qmax_i", "rhomax_i", "L_j", "v_j", "w_j", "qmax_j", "rhomax_j", "integral", "pi"]
        writer.writerow(header)
        with multiprocessing.Pool(processes=os.cpu_count()) as pool:
            # pool = multiprocessing.Pool() #generate processes equal to the number of cores
            chunksize = ceil(vaex.dataframe.DataFrame.count(dv) / os.cpu_count())
            #result = p_map(mm.simulation, dv.values)

            for result in pool.map(mm.simulation, dv.values, chunksize=chunksize):
                writer.writerow(result)

    elapsed = time.time() - t
    print("Elapsed time: " + str(elapsed))
    # aa = np.array(result[0])
    # print(aa)
    # dv['integral'] = np.array(result[0])
    # dv['pi'] = np.array(result[1])

     #dv.export_csv('C:/Users/User/Desktop/cazzillo/output_sim1.csv', sep=';')
    #dv.export_hdf5('C:/Users/User/Desktop/cazzillo/output_sim1.hdf5')


