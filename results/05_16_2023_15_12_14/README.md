dt = 10
window_length = 270
freq = 30

congestion_start = 2400
congestion_end = 4000

Note: Results from 05_12_2023 onwards are with the true epigraph reformulation
Note: Results from 05_16_2023 onwards are with better results saving, also saves supply and demand

Note: Longer and larger congestion period
phi_day = 1.2 * gen.get_flow(t_0=day_start, t_f=day_end, perturb=False)

class ControlParameters:
    self.discount = 0.05
    self.epi_st = 0.1
    self.a_rho = 0.01
    self.a_queue = 0.0009
    self.a_epi = 0.003
    self.a_cont = 0.70