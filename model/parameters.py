from dataclasses import dataclass
from typing import Dict
import numpy as np
import pandas as pd


# TODO: Add assert statements to check that data is well formatted
# TODO: Add def __post_init__(self)
# TODO: Make sure that @dataclass is the right option. More of a python question than anything
@dataclass
class HighwayParameters:
    """
    Class to represent all parameters for the highway in the CTM-s model
    """

    id: np.ndarray          # Cell ID [0,1,...,N-1]
    l: np.ndarray           # Cell Length [km]
    v: np.ndarray           # Free-flow Speed [km/hr]
    w: np.ndarray           # Congestion Wave Speed [km/hr]
    q_max: np.ndarray       # Maximum Cell Capacity [veh/hr]
    rho_max: np.ndarray     # Maximum Jam Density [veh/km]
    p_ms: np.ndarray        # Mainstream priority of each cell [0,1]
    dt: float               # dt [hrs]

    def __init__(self, loc: str):
        parameters = pd.read_csv(loc, sep=';').to_numpy()
        self.id = parameters[:, 0].astype(int)
        self.l = parameters[:, 1]
        self.v = parameters[:, 2]
        self.w = parameters[:, 3]
        self.q_max = parameters[:, 4]
        self.rho_max = parameters[:, 5]
        self.p_ms = np.ones(len(self.id))
        self.dt = 10/3600  # Default value (10 seconds)

    def __len__(self):
        return len(self.id)


@dataclass
class OnRampParameters:
    """
    Class to represent all parameters for on-ramps in the CTM-s model
    """

    id: np.ndarray          # On-ramp ID [0,1,...,N-1]
    j: np.ndarray           # On-ramp Access Cell
    r_r_max: np.ndarray     # Maximum On-ramp Flow [veh/hr]
    p_r: np.ndarray         # Priority of On-ramp [0,1]
    d_r: np.ndarray         # External Flow Signal [veh/hr]

    def __init__(self, loc: str, phi_loc: str):
        parameters = pd.read_csv(loc, sep=';').to_numpy()
        self.id = parameters[:, 0].astype(int)
        self.j = parameters[:, 1].astype(int)
        self.r_r_max = parameters[:, 2]
        self.p_r = parameters[:, 3]

        if len(self.id) == 0:
            self.d_r = np.empty(0)
        else:
            self.d_r = pd.read_csv(phi_loc, header=None, sep=';').to_numpy()

    def __len__(self):
        return len(self.id)


@dataclass
class OffRampParameters:
    """
    Class to represent all parameters for off-ramps in the CTM-s model
    """

    id: np.ndarray          # Off-ramp ID [0,1,...,N-1]
    i: np.ndarray           # Off-ramp Exit Cell
    beta_r: np.ndarray      # Split Ratio [0,1]

    def __init__(self, loc: str):
        parameters = pd.read_csv(loc, sep=';').to_numpy()
        self.id = parameters[:, 0].astype(int)
        self.i = parameters[:, 1].astype(int)
        self.beta_r = parameters[:, 2]

    def __len__(self):
        return len(self.id)


@dataclass
class StationParameters:
    """
    Class to represent all parameters for stations in the CTM-s model
    """

    id: np.ndarray          # Station ID [0,1,...,N-1]
    i: np.ndarray           # Station Access Cell
    j: np.ndarray           # Station Exit Cell
    j_r: np.ndarray         # Unique Station Exit Cells
    r_s_max: np.ndarray     # Maximum Service Station Exit Flow [veh/hr]
    delta: np.ndarray       # Average Time Spent at Service Station [Time increments]
    beta_s: np.ndarray      # Split Ratio (Entering Service Station) [0,1]
    p: np.ndarray           # Priority (Exiting Service Station) [0,1]
    j_to_p: Dict            # Map: j -> p; cell j -> priority of all stations exiting at cell j

    def __init__(self, loc: str):
        station = pd.read_csv(loc, sep=';').to_numpy()
        self.id = station[:, 0].astype(int)
        self.i = station[:, 1].astype(int)
        self.j = station[:, 2].astype(int)
        self.j_r = np.unique(self.j).astype(int)
        self.r_s_max = station[:, 3]
        self.delta = station[:, 4]
        self.beta_s = station[:, 5]  # Default value
        self.p = station[:, 6]

        self.build_j_to_p()

    def build_j_to_p(self):
        """
        Build map from station outlet cell to total priority
        """

        self.j_to_p = {}

        for j in self.j_r:
            self.j_to_p[j] = 0

        for i_, j in enumerate(self.j):
            self.j_to_p[j] += self.p[i_]

    def __len__(self):
        return len(self.id)


class CTMsParameters:
    """
    Class to represent all parameters in the CTM-s model
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

        self.update_mainstream_priorities()

    def update_mainstream_priorities(self):
        """
        Calculate mainstream_priorities by considering both stations ond on-ramps
        """

        for j_, j in enumerate(self.onramps.j):
            self.highway.p_ms[j] -= self.onramps.p_r[j_]

        for j_, j in enumerate(self.stations.j):
            self.highway.p_ms[j] -= self.stations.p[j_]
