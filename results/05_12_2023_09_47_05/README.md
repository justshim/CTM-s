dt = 10
window_length = 270
freq = 60

congestion_start = 2400
congestion_end = 4000

Note: Results from 05_12_2023 and moving forward are with the true epigraph reformulation

Note: Longer and larger congestion period
phi_day = 1.2 * gen.get_flow(t_0=day_start, t_f=day_end, perturb=False)

class ControlParameters:
    alpha_epi = 0.003
    alpha_tts = 0.70

a_epigraph = -np.ones(self.num_variables)
a_epigraph[16::self.m] = -0.7
discount = np.kron(np.exp(-0.05 * np.arange(0, self.k_l)), np.ones(self.m))
a_epigraph[:-1] = np.multiply(a_epigraph[:-1], discount)
a_epigraph = np.atleast_2d(a_epigraph)
b_epigraph = np.zeros((1, 1))

self.add_constraints(a_epigraph, b_epigraph)

def define_cost_function(self):
    tts_1 = 1 * np.sum(self.rho_matrix, axis=0)
    tts_2 = 0.3 * np.sum(self.e_matrix, axis=0)
    tts = self.dt * (tts_1 + tts_2)  # Units: [veh hr / km]
    discount = np.kron(np.exp(-0.05 * np.arange(0, self.k_l)), np.ones(self.m))
    tts[:-1] = np.multiply(tts[:-1], discount)
    epi = np.zeros(self.num_variables)
    epi[-1] = 1
    self.c = (self.params_c.alpha_epi * epi) + (self.params_c.alpha_tts * tts)