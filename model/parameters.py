import numpy as np
import pandas as pd


class HighwayParameters:
    ID: np.ndarray          # Cell ID [0,1,...,N-1]
    L: np.ndarray           # Cell Length [km]
    v: np.ndarray           # Free-flow Speed [km/hr]
    w: np.ndarray           # Congestion Wave Speed [km/hr]
    q_max: np.ndarray       # Maximum Cell Capacity [veh/hr]
    rho_max: np.ndarray     # Maximum Jam Density [veh/km]
    dt: np.ndarray          # dt [] TODO: Units?

    def __init__(self, loc: str):
        parameters = pd.read_csv(loc, sep=';').to_numpy()
        self.ID = parameters[:, 0]
        self.L = parameters[:, 1]
        self.v = parameters[:, 2]
        self.w = parameters[:, 3]
        self.q_max = parameters[:, 4]
        self.rho_max = parameters[:, 5]
        self.dt = parameters[:, 6]

    def __len__(self):
        return len(self.ID)


class OnRampParameters:
    ID: np.ndarray          # On-ramp ID [0,1,...,N-1]
    j: np.ndarray           # On-ramp Access Cell
    r_r_max: np.ndarray     # Maximum On-ramp Flow [veh/hr]
    p_r: np.ndarray         # Priority of On-ramp [0,1]
    d_r: np.ndarray         # External Flow Signal [veh/hr]

    def __init__(self, loc: str, phi_loc: str):
        parameters = pd.read_csv(loc, sep=';').to_numpy()
        self.ID = parameters[:, 0]
        self.j = parameters[:, 1]
        self.r_r_max = parameters[:, 2]
        self.p_r = parameters[:, 3]
        self.d_r = pd.read_csv(phi_loc, sep=';').to_numpy()

    def __len__(self):
        return len(self.ID)


class OffRampParameters:
    ID: np.ndarray          # Off-ramp ID [0,1,...,N-1]
    i: np.ndarray           # Off-ramp Exit Cell
    beta_r: np.ndarray      # Split Ratio [0,1]

    def __init__(self, loc: str):
        parameters = pd.read_csv(loc, sep=';').to_numpy()
        self.ID = parameters[:, 0]
        self.i = parameters[:, 1]
        self.beta_r = parameters[:, 2]

    def __len__(self):
        return len(self.ID)


class StationParameters:
    ID: np.ndarray          # Station ID [0,1,...,N-1]
    i: np.ndarray           # Station Access Cell
    j: np.ndarray           # Station Exit Cell
    r_s_max: np.ndarray     # Maximum Service Station Exit Flow [veh/hr]
    delta: np.ndarray       # Average Time Spent at Service Station []  # TODO: Units?
    beta_s: np.ndarray      # Split Ratio (Entering Service Station) [0,1]
    p: np.ndarray           # Priority (Exiting Service Station) [0,1]

    def __init__(self, loc: str):
        station = pd.read_csv(loc, sep=';').to_numpy()
        self.ID = station[:, 0]
        self.i = station[:, 1]
        self.j = station[:, 2]
        self.r_s_max = station[:, 3]
        self.delta = station[:, 4]
        self.beta_s = station[:, 5]
        self.p = station[:, 6]

    def __len__(self):
        return len(self.ID)


class CTMsParameters:
    """
    Basic Class to represent all CTMs Parameters
    """
    highway: HighwayParameters
    onramps: OnRampParameters
    offramps: OffRampParameters
    stations: StationParameters

    def __init__(self, hi_loc: str, onr_loc: str, offr_loc: str, st_loc: str, phi_onr_loc: str):
        self.highway = HighwayParameters(hi_loc)
        self.onramps = OnRampParameters(onr_loc, phi_onr_loc)
        self.offramps = OffRampParameters(offr_loc)
        self.stations = StationParameters(st_loc)
