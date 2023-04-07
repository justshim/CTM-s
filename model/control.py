import numpy as np


# TODO: Add Control Parameter Class
# class TrafficControlParameters:
#     raise NotImplementedError

class TrafficControlInput:
    """
    TODO: ...
    """
    # TODO: For all k in time interval
    r_s_c: np.ndarray           # Service Station Off-Ramp Metering [veh/hr]
    beta_s: np.ndarray          # Service Station On-Ramp Split Ratio [0,1] # TODO: k_interval x num_stations
