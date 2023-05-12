dt = 10
window_length = 270
freq = 60

congestion_start = 2400
congestion_end = 4000

def define_cost_function(self):
    tts_1 = 1 * np.sum(self.rho_matrix, axis=0)
    tts_2 = 0.3 * np.sum(self.e_matrix, axis=0)
    tts = self.dt * (tts_1 + tts_2)  # Units: [veh hr / km]
    ttd = self.dt * np.ones(self.num_variables)
    self.c = (0.7 * tts) - (1 * ttd)  # TODO: self.eta = 1 
    self.c[::self.m] = -1 
    discount = np.kron(np.exp(-0.05 * np.arange(0, self.k_l)), np.ones(self.m))
    self.c = np.multiply(self.c, discount)

Note: Longer and larger congestion period
phi_day = 1.2 * gen.get_flow(t_0=day_start, t_f=day_end, perturb=False)

Note: Less prioritization on maximizing the service station outflow
TODO: By doing the proper epigraph reformulation, this would be easier...
ttd_1 = np.ones(self.num_variables)
ttd_1[16::self.m] *= 0.3  
ttd = self.dt * ttd_1