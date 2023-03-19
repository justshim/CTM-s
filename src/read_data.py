from model.parameters import *


def read_parameters(hi_loc: str, onr_loc: str, offr_loc: str, st_loc: str, phi_onr_loc: str) -> CTMsParameters:
    return CTMsParameters(hi_loc, onr_loc, offr_loc, st_loc, phi_onr_loc)


def read_phi(loc: str) -> np.ndarray:
    return pd.read_csv(loc, header=None).to_numpy().flatten()