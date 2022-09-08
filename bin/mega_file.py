import itertools
import random
import numpy as np
import pandas as pd
import vaex

from model.optimization import main_mega as mm


if __name__ == '__main__':
    dv = vaex.open('C:/Users/User/Desktop/cazzillo/output_data_random.hdf5')
    path_phi = ''

    row = dv.take(1)
    print(row)
    # for idx in range(vaex.dataframe.DataFrame.count(dv)):
    #     row = dv.take(idx)
    #     mm.simulation(dv, path_phi)


