from typing import List
from enum import IntEnum


class CongState(IntEnum):
	"""
	Four congestion states as described in
	"A Novel Control-Oriented Cell Transmission Model Including Service Stations on Highways"
	Cenedese et al. 2022
	"""

	FREEFLOW = 0
	CONG_MS = 1
	CONG_ST = 2
	CONG_ALL = 3


class Cell:
	"""
	Class modeling the cells on a highway stretch in the CTM-s model
	"""

	id: int								# Cell ID
	l: float							# Length of Cell [km]
	v_free: float						# Freeflow Velocity [km/hr]
	w: float							# Congestion Wave Velocity [km/hr]
	p_ms: float							# Mainstream Priority [0,1]
	q_max: float						# Maximum Supported Vehicle Flow [veh/hr]
	rho_max: float						# Maximum Supported Vehicle Density [veh/km]

	rho: List[float]					# Cell Vehicle Density [veh/km]
	phi: List[float]					# Cell Incoming Mainstream Vehicle Flow [veh/hr]
	v: List[float]						# Cell Velocity [km/hr]
	phi_minus: List[float]				# Total Exit Vehicle Flow [veh/hr]
	phi_plus: List[float]				# Total Entering Vehicle Flow [veh/hr]

	dem: List[float]					# Cell Demand [veh/hr]
	sup: List[float]					# Cell Supply [veh/hr]
	cong: List[CongState]				# Cell Congestion State

	k: int								# Time step

	def __init__(self, id_cell, length, v_free, w, q_max, rho_max, p_ms):
		self.id = id_cell
		self.l = length
		self.v_free = v_free
		self.w = w
		self.p_ms = p_ms
		self.q_max = q_max
		self.rho_max = rho_max

		self.rho = [0]
		self.phi = []
		self.v = []
		self.phi_minus = []
		self.phi_plus = []

		self.dem = []
		self.sup = [0]
		self.cong = []

		self.k = 0

	def to_string(self):
		"""
		Utility method to print some information about the cell
		"""

		print("Cell ID: " + str(self.id))
		print("Length: " + str(self.l))
		print("v_free: " + str(self.v_free))
		print("w: " + str(self.w))
		print("rho_max: " + str(self.rho_max))
		print("q_max: " + str(self.q_max))
		print()

	def compute_velocity(self):
		"""
		Compute cell velocity

		Note: CTM-s is only a 1st order model, so v <= v_free is not respected for all times
		"""

		if self.rho[self.k] == 0:
			self.v.append(self.v_free)
		else:
			phi_avg = self.phi_plus[-1]
			self.v.append(phi_avg/self.rho[self.k])

	def compute_phi(self, d_prev: float, total_ds: float):
		"""
		Computation of the flow entering this cell from the previous one at time instant k
		Varies according to the congestion state of the cell, so we update it first
		"""

		self.update_cong_state(d_prev, total_ds)

		cong_state = self.cong[self.k]
		sup = self.sup[self.k]

		if cong_state == CongState.FREEFLOW:
			self.phi.append(d_prev)
		elif cong_state == CongState.CONG_MS:
			self.phi.append(sup - total_ds)
		elif cong_state == CongState.CONG_ST:
			self.phi.append(d_prev)
		elif cong_state == CongState.CONG_ALL:
			self.phi.append(sup * self.p_ms)

	def compute_phi_plus(self, rs_total: float):
		"""
		Computation of the total flow entering this cell at time instant k
		"""

		self.phi_plus.append(self.phi[self.k] + rs_total)
		
	def compute_phi_minus(self, ss: float, next_phi: float):
		"""
		Computation of the total flow exiting this cell at time instant k
		"""

		self.phi_minus.append(next_phi + ss)

	def compute_rho(self, dt: float):
		"""
		Computation of the traffic density of this cell at time instant k + 1
		"""

		self.rho.append(self.rho[self.k] + (dt/self.l) * (self.phi_plus[self.k] - self.phi_minus[self.k]))

	def compute_demand(self, total_beta: float):
		"""
		Computation of the demand of this cell at time instant k
		"""

		self.dem.append(min([(1 - total_beta) * self.v_free * self.rho[self.k], self.q_max]))
		
	def compute_supply(self):
		"""
		Computation of the supply of this cell at time instant k
		"""

		self.sup.append(min([self.w * (self.rho_max - self.rho[self.k]), self.q_max]))

	def update_cong_state(self, d_prev: float, total_ds: float):
		"""
		Computation of the congestion state of this cell at time instant k
		"""

		sup = self.sup[self.k]

		if d_prev + total_ds <= sup:
			self.cong.append(CongState.FREEFLOW)
		
		elif (d_prev > self.p_ms * sup) and (total_ds <= (1 - self.p_ms) * sup):
			self.cong.append(CongState.CONG_MS)

		elif (d_prev <= self.p_ms * sup) and (total_ds > (1 - self.p_ms) * sup):
			self.cong.append(CongState.CONG_ST)

		elif (d_prev > self.p_ms * sup) and (total_ds > (1 - self.p_ms) * sup):
			self.cong.append(CongState.CONG_ALL)
		else:
			print("cell: " + str(self.id) + " Congestion Error")

	def update_k(self, kappa: int):
		"""
		Each iteration starts with the update of the time instant
		"""

		self.k = kappa

	def reset(self):
		"""
		Reset cell to zero initial condition
		"""

		self.v = []
		self.phi = []
		self.phi_minus = []
		self.phi_plus = []
		self.rho = [0]
		self.dem = []
		self.sup = [0]
		self.cong = []
		self.k = 0
