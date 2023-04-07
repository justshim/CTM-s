from model.parameters import *


def read_parameters(hi_loc: str, onr_loc: str, offr_loc: str, st_loc: str, phi_onr_loc: str) -> CTMsParameters:
    """
    Helper function to read in CTM-s Parameters
    """

    return CTMsParameters(hi_loc, onr_loc, offr_loc, st_loc, phi_onr_loc)


def read_phi(loc: str, t_0: int, t_f: int, ds=1) -> np.ndarray:
    """
    Read in real-flow data from t_0 to t_f
    Flow data is measured every 10 seconds over the course of 24 hours
    ds = down-sampling rate e.g. ds = 6 => flow every minute
    """

    phi = pd.read_csv(loc, header=None).to_numpy().flatten()
    return phi[t_0:t_f:int(ds)]
