import itertools
import random
import numpy as np
import pandas as pd
import vaex

from model.optimization import main_mega as mm


if __name__ == '__main__':
    dv = vaex.open('C:/Users/User/Desktop/cazzillo/output_data_random.hdf5')
    path_phi = 'C:/Users/User/Documents/MATLAB/CTMs-identification/fnc/extracted_data/CTM_param_out_nice.xls'
    #print(dv)
    # row = dv.take([1])
    # a = row.values.flatten()
    # print(row)
    results = pd.DataFrame({"integral": np.empty(vaex.dataframe.DataFrame.count(dv)), "pi": np.empty(vaex.dataframe.DataFrame.count(dv))})
    #print(results)
    for idx in range(2):
        row = dv.take([idx])
        #print(row.values.flatten())
        results.loc[idx] = mm.simulation(row.values.flatten(), path_phi)

    #print(results['integral'].astype(str))
    dv['integral'] = results['integral'].astype(str)
    dv['pi'] = results['pi'].astype(str)
    print(dv)



