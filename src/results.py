from typing import List

import numpy as np
from matplotlib import pyplot as plt

# class MyTrafficResults:
#     def __init__(self, day_length: int):
#         self.day_length = day_length


def plot_comparison(y_0: np.ndarray, y_1: np.ndarray, i: List, k_0: int, k_f: int, loc: str):
    fig = plt.figure()
    t = np.arange(k_0, k_f)
    plt.plot(t, y_0[:, i], '-k')
    plt.plot(t, y_1[:, i], '--r')
    plt.xlabel('Time')
    fig.savefig(loc + 'comp/cell.png', dpi=300)  # TODO: Naming scheme for files...
