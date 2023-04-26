from typing import List
import numpy as np

from model.cell import Cell
from model.station import Station
from model.on_ramp import OnRamp
from model.off_ramp import OffRamp


class Stretch:
	"""
	Class to maintain and simulate CTM-s model
	"""

	cells: List[Cell]			# Cell in highway stretch
	on_ramps: List[OnRamp]		# On-ramps in highway stretch
	off_ramps: List[OffRamp]  	# Off-ramps in highway stretch
	stations: List[Station]		# Stations in highway stretch

	n_cells: int				# Number of Cells
	n_on_ramps: int				# Number of On-ramps
	n_off_ramps: int			# Number of Off-ramps
	n_stations: int				# Number of Stations
	n_stations_exit: int		# Number of Station Exit Flow Points

	phi_0: np.ndarray			# Incoming flow of first cell [veh/hr]
	dt: float					# Time increment [hrs]
	k: int						# Time step

	y_rho: np.ndarray			# State: Density of Cells
	y_l: np.ndarray				# State: Number of Service Station Users
	y_e: np.ndarray				# State: Length of Service Station Queue

	u_phi: np.ndarray			# Input: Flows between cells
	u_rs: np.ndarray			# Input: Flow exiting service station

	s_s: np.ndarray				# Flow entering service station (for LP initialization)
	beta_total: np.ndarray		# Aggregate split ratio at each cell

	def __init__(self, dt: float):
		self.cells = []
		self.on_ramps = []
		self.off_ramps = []
		self.stations = []

		self.n_cells = 0
		self.n_on_ramps = 0
		self.n_off_ramps = 0
		self.n_stations = 0
		self.n_stations_exit = 0

		self.phi_0 = np.empty(0)
		self.dt = dt
		self.k = 0

		self.y_rho = np.empty(0)
		self.y_l = np.empty(0)
		self.y_e = np.empty(0)

		self.u_phi = np.empty(0)
		self.u_rs = np.empty(0)
		self.s_s = np.empty(0)
		self.beta_total = np.empty(0)

	def to_string(self):
		"""
		Utility method to print some information about the highway stretch
		"""

		for cell in self.cells:
			cell.to_string()

		for on_ramp in self.on_ramps:
			on_ramp.to_string()

		for off_ramp in self.off_ramps:
			off_ramp.to_string()

		for station in self.stations:
			station.to_string()

		print("Number of cells: " + str(self.n_cells))
		print("Number of on-ramps: " + str(self.n_on_ramps))
		print("Number of off-ramps: " + str(self.n_off_ramps))
		print("Number of stations: " + str(self.n_stations))

		print("Time Length: " + str(self.dt))
		print()
	
	def create_cell(self, length, v_free, w, q_max, rho_max, p):
		"""
		Method to create an instance of the object Cell, and add it to this stretch
		"""

		cell = Cell(self.n_cells, length, v_free, w, q_max, rho_max, p)
		self.cells.append(cell)
		self.n_cells += 1

	def create_on_ramp(self, d_r, r_r_max, j, p_r):
		"""
		Method to create an instance of the object Station, and add it to this stretch
		"""

		on_ramp = OnRamp(self.n_on_ramps, d_r, r_r_max, j, p_r)
		self.on_ramps.append(on_ramp)
		self.n_on_ramps += 1

	def create_off_ramp(self, i, beta_r):
		"""
		Method to create an instance of the object Station, and add it to this stretch
		"""

		off_ramp = OffRamp(self.n_off_ramps, i, beta_r)
		self.off_ramps.append(off_ramp)
		self.n_off_ramps += 1

	def create_station(self, r_s_max, i, j, delta, beta_s, p):
		"""
		Method to create an instance of the object Station, and add it to this stretch
		"""

		if j not in [s.j for s in self.stations]:
			self.n_stations_exit += 1

		station = Station(self.n_stations, r_s_max, i, j, delta, beta_s, p)
		self.stations.append(station)
		self.n_stations += 1

	def update(self, k):
		"""
		Main method of the class
		At each time instant k, updates all the parameters of the cells and service stations on this stretch
		"""

		# Update time instant
		self.k = k

		# Preliminary updates
		# First update time instant for all cells with current k
		for i in range(len(self.cells)):
			self.cells[i].update_k(self.k)

		# Same for stations, plus computation of some preliminary values
		for s in range(len(self.stations)):
			self.stations[s].update_k(self.k)
			self.stations[s].compute_demand(self.dt)

		# Same for on-ramps, plus computation of some preliminary values
		for r_on in range(len(self.on_ramps)):
			self.on_ramps[r_on].update_k(self.k)
			self.on_ramps[r_on].computeDrBig(self.dt)
		
		# First batch of cell value updates, with special case for first cell
		for i in range(len(self.cells)):
			# Compute demand from cell i for cell i+1 (phi_i+1)
			self.beta_total[self.k, i] = self.compute_total_beta(i)
			self.cells[i].compute_demand(self.beta_total[self.k, i])

			# Compute demand from cell i-1 for cell i (phi_i)
			total_ds = self.compute_total_ds(i)
			d_prev = self.compute_prev_demand(i)  # Cell 0: prev. demand as input (first_phi)
			self.cells[i].compute_phi(d_prev, total_ds)

		# Second batch of cell value updates, with special case for last cell
		for i in range(len(self.cells)):  # N = Number of Cells
			# Cell N-1: phi_N is the flow leaving cell N-1
			# phi_N = Demand from cell N-2 (phi_N-1)
			next_phi = self.compute_next_phi(i)

			total_rs_station = self.compute_r_station(i)
			total_ss_station = self.compute_s_station(i, self.beta_total[self.k, i], next_phi)

			total_sr_ramp = self.compute_s_ramp(i, next_phi)
			total_rr_ramp = self.compute_r_ramp(i)

			total_rs = total_rr_ramp + total_rs_station
			total_ss = total_sr_ramp + total_ss_station

			self.cells[i].compute_supply()
			self.cells[i].compute_phi_minus(total_ss, next_phi)
			self.cells[i].compute_phi_plus(total_rs)

		# Final updates
		# All cells have their densities updated
		for i in range(len(self.cells)):
			self.cells[i].compute_rho(self.dt)

		# All stations have their l and e updated
		for s in range(len(self.stations)):
			self.stations[s].compute_queue_length(self.dt)
			self.stations[s].compute_num_vehicles(self.dt)

		# All ramps have their l updated
		for r_on in range(len(self.on_ramps)):
			self.on_ramps[r_on].computeL(self.dt)

		# Compute cell velocities
		for i in range(len(self.cells)):
			self.cells[i].compute_velocity()

		# Track Variables
		for i, c in enumerate(self.cells):
			self.y_rho[self.k + 1, i] = c.rho[-1]
			self.u_phi[self.k, i] = c.phi[-1]

		for i, s in enumerate(self.stations):
			self.y_l[self.k + 1, i] = s.l[-1]
			self.y_e[self.k + 1, i] = s.e[-1]
			self.u_rs[self.k, i] = s.r_s[-1]
			self.s_s[self.k, i] = s.s_s[-1]

	def compute_total_beta(self, i: int) -> float:
		"""
		Compute total beta from service-stations and off-ramps, needed to compute D_i
		"""

		total_beta = 0

		# Check for service-stations with access point at cell i
		for s in range(len(self.stations)):
			if self.stations[s].i == i:
				total_beta += self.stations[s].beta_s

		# Check for off-ramps exiting at cell i
		for r_off in range(len(self.off_ramps)):
			if self.off_ramps[r_off].i == i:
				total_beta += self.off_ramps[r_off].beta_r

		return total_beta

	def compute_total_ds(self, i: int) -> float:
		"""
		Compute total D_s from service-stations and on-ramps, needed to compute phi_i
		"""

		total_ds = 0

		# Check for service-stations with exit point at cell i
		for s in range(len(self.stations)):
			if self.stations[s].j == i:
				total_ds += self.stations[s].demand

		# Check for on-ramps entering at cell i
		for r_on in range(len(self.on_ramps)):
			if self.on_ramps[r_on].j == i:
				total_ds += self.on_ramps[r_on].d_r_big

		return total_ds
	
	def compute_r_station(self, i: int) -> float:
		"""
		For each cell, check if any stations merge into it, and compute their r_s
		"""

		total_rs = 0
		for s in range(len(self.stations)):
			if self.stations[s].j == i:
				if self.cells[i].cong_state == 0 or self.cells[i].cong_state == 1:
					self.stations[s].compute_rs(0, self.cells[i].cong_state)
					
				elif self.cells[i].cong_state == 2:
					self.iterative_procedure(i, self.cells[i].cong_state)
					
				elif self.cells[i].cong_state == 3:
					self.iterative_procedure(i, self.cells[i].cong_state)
					
				total_rs += self.stations[s].r_s[-1]

		return total_rs

	def iterative_procedure(self, i: int, cong_state: int):
		"""
		Method called during the update procedure and used to assign r_s to all stations merging into the same cell
		in case of congestions of type 2 and 3
		"""

		# initialization of support variables
		demands = []
		prev_demand = self.cells[i - 1].demand
		supply = self.cells[i].supply
		good = [0]  # list to contain "good" demands, i.e. the ones that do not saturate the flow
		sum_d_good = 0
		sum_p = 0

		for s in self.stations:
			if s.j == i:
				demands.append(s)

		if cong_state == 2:
			supply_res = supply - prev_demand

		elif cong_state == 3:
			supply_res = (1 - self.cells[i].p_ms) * supply

		bad = demands  # list to contain "bad" demands, i.e. the ones that do saturate the flow

		# Recursively compute "bad" (E_cal_overline in the paper) and "good" (E_cal_underline in the paper)
		while len(good) != 0:
			good.clear()
			if cong_state == 2:
				for d in demands:
					if d.demand <= (supply_res - sum_d_good) / len(bad):
						bad.remove(d)
						good.append(d)

						# update RS for "good"
						for station in self.stations:
							if d.id_station == int(station.id_station):
								station.compute_rs(d.demand, cong_state)

						sum_d_good = sum_d_good + d.demand

						supply_res = supply_res - d.demand
			elif cong_state == 3:
				for d in demands:
					if d.demand <= (((1 - self.cells[i].p_ms) * supply_res) - sum_d_good) / len(bad):
						bad.remove(d)
						good.append(d)

						# update RS for "good"
						for station in self.stations:
							if d.id_station == int(station.id_station):
								station.compute_rs(d.demand, cong_state)

						sum_d_good = sum_d_good + d.demand

						supply_res = (1 - self.cells[i].p_ms) * supply_res - d.demand  ## VERIFICARE FORMULA

		# Compute sum of priorities for all involved stations
		for b in bad:
			sum_p = sum_p + b.p

		# Update RS for "bad"
		for b in bad:
			for station in self.stations:
				if b.id_station == int(station.id_station):
					station.compute_rs((b.p / sum_p) * supply_res, cong_state)

	def compute_s_station(self, i: int, beta_total: float, next_phi: float) -> float:
		"""
		Compute total flow exiting cell i to service stations
		"""

		total_ss = 0
		for s in range(len(self.stations)):
			if self.stations[s].i == i:
				self.stations[s].compute_ss(beta_total, next_phi)
				total_ss += self.stations[s].s_s[self.k]

		return total_ss

	def compute_r_ramp(self, i: int) -> float:
		"""
		Compute total flow entering cell i via on-ramp
		"""

		total_rs = 0
		for r_on in range(len(self.on_ramps)):

			if self.on_ramps[r_on].j == i:
					
				if self.cells[i].cong_state == 0 or self.cells[i].cong_state == 1:
					self.on_ramps[r_on].computeRr(0, self.cells[i].cong_state)
				
				elif self.cells[i].cong_state == 2:
					rr = self.cells[i].supply - self.cells[i - 1].demand
					self.on_ramps[r_on].computeRr(rr, self.cells[i].cong_state)
				
				elif self.cells[i].cong_state == 3:
					rr = self.cells[i].supply * self.on_ramps[r_on].p_r
					self.on_ramps[r_on].computeRr(rr, self.cells[i].cong_state)
				
				total_rs += self.on_ramps[r_on].r_r

		return total_rs

	def compute_s_ramp(self, i: int, next_phi: float) -> float:
		"""
		Compute total flow exiting cell i via off-ramps
		"""

		total_ss = 0
		for r_off in range(len(self.off_ramps)):

			if self.off_ramps[r_off].i == i:
				self.off_ramps[r_off].computeSr(next_phi)
				total_ss += self.off_ramps[r_off].s_r

		return total_ss

	def compute_prev_demand(self, i: int) -> float:
		"""
		Compute demand of cell i-1, needed to determine phi_i
		Note: First cell does not have a "previous" cell, hence phi_(i-1) is given as input
		"""

		if i != 0:
			prev_demand = self.cells[i-1].demand
		else:
			prev_demand = self.phi_0[self.k]
			
		return prev_demand

	def compute_next_phi(self, i: int) -> float:
		"""
		Compute flow leaving cell i, phi_(i+1)
		Note: Last cell does not have a "next" cell, so flow exiting last cell set to demand of penultimate cell
		Result: Flow entering last cell equals flow exiting last cell, last cell density = 0 for all time
		"""

		if i < len(self.cells) - 1:
			next_phi = self.cells[i+1].phi[-1]
		else:
			# next_phi = self.last_phi[self.k]
			next_phi = self.cells[i-1].demand
			self.u_phi[self.k, -1] = next_phi

		return next_phi

	def simulate_init(self, phi_0: np.ndarray):
		"""
		Initialize arrays to track the values of states and inputs during simulation
		"""

		self.phi_0 = phi_0
		day_length = len(phi_0)
		self.y_rho = np.zeros((day_length + 1, self.n_cells))
		self.y_l = np.zeros((day_length + 1, self.n_stations))
		self.y_e = np.zeros((day_length + 1, self.n_stations))

		self.u_phi = np.zeros((day_length, self.n_cells + 1))
		self.u_rs = np.zeros((day_length, self.n_stations_exit))
		self.s_s = np.zeros((day_length, self.n_stations))
		self.beta_total = np.zeros((day_length, self.n_cells))

		# Default control value
		for s in self.stations:
			s.r_s_c = s.r_s_max * np.ones((day_length, 1))

	def get_state(self, k: int) -> np.ndarray:
		"""
		Get state, as defined in the LP, at time instant k
		"""

		state = np.zeros(self.n_cells + 2*self.n_stations)
		state[0:self.n_cells] = self.y_rho[k, :]
		state[self.n_cells::2] = self.y_l[k, :]
		state[self.n_cells+1::2] = self.y_e[k, :]

		return state.reshape((-1, 1))

	def get_station_inflow(self, k: int) -> np.ndarray:
		"""
		Get history of flow into service station until time instant k
		"""

		return self.s_s[:k, :]

	def update_control(self, k_0: int, r_s_c: np.ndarray):
		"""
		Update ramp-meter control for each service station on the highway from time k_0
		"""

		k_l = len(r_s_c)
		for i in range(self.n_stations):
			self.stations[i].r_s_c[k_0:k_0 + k_l, i] = r_s_c[:, i]

	def reset(self):
		"""
		Reset all objects in the stretch
		"""

		self.k = 0

		for cell in self.cells:
			cell.reset()

		for on_ramp in self.on_ramps:
			on_ramp.reset()

		for off_ramp in self.off_ramps:
			off_ramp.reset()

		for station in self.stations:
			station.reset()
