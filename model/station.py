from typing import List
import numpy as np

from model.cell import CongState


class Station:
	"""
	Class modeling the service stations on a highway stretch in the CTM-s model
	"""

	id: int 					# Service Station ID
	i: int						# Service Station Inlet Cell ID
	j: int						# Service Station Outlet Cell ID
	delta: float				# Time Spent at Service Station [Time Steps]
	beta_s: float				# Service Station Inflow Split Ratio [0,1]
	p: float					# Service Station Outflow Priority [0,1]
	r_s_max: float				# Maximum Supported Outflow [veh/hr]
	e_max: float  				# TODO: Maximum Queue Length

	e: List[float]				# Service Station Queue Length [veh]
	l: List[float]				# Number of Service Station Users [veh]

	s_s: List[float]			# Service Station On-ramp Flow [veh/hr]
	s_e: List[float]			# Station to Queue Flow [veh/hr]
	r_s: List[float]			# Service Station Off-ramp Flow [veh/hr]

	dem: List[float]			# Station Off-ramp Demand [veh/hr]

	r_s_c: np.ndarray			# Service Station Off-ramp Control Flow [veh/hr]

	k: int						# Time step

	def __init__(self, id_station: int, i: int, j: int, delta: float, beta_s: float, p: float, r_s_max: float):
		self.id = id_station
		self.i = i
		self.j = j
		self.delta = delta
		self.beta_s = beta_s
		self.p = p
		self.r_s_max = r_s_max
		self.e_max = 100  # TODO: Hardcoded

		self.e = [0]
		self.l = [0]

		self.s_s = []
		self.s_e = []
		self.r_s = []

		self.dem = []

		self.r_s_c = np.empty(0)  # Default value

		self.k = 0

	def to_string(self):
		"""
		Utility method to print some information about the service station
		"""

		print("Station ID: " + str(self.id))
		print("From cell " + str(self.i))
		print("To cell " + str(self.j))
		print("Time delay: " + str(self.delta))
		print("Split ratio: " + str(self.beta_s))
		print("Number of vehicles: " + str(self.l))
		print()

	def compute_ss(self, beta_total: float, next_phi: float):
		"""
		Computation of the flow leaving the mainstream to enter this service station at time instant k
		"""

		self.s_s.append((self.beta_s / (1 - beta_total)) * next_phi)
		
	def compute_rs(self, rs: float, cong_state: CongState):
		"""
		Computation of the flow merging into the mainstream from this service station at time instant k
		"""

		if cong_state == CongState.FREEFLOW or cong_state == CongState.CONG_MS:
			self.r_s.append(self.dem[self.k])
		elif cong_state == CongState.CONG_ST or cong_state == CongState.CONG_ALL:
			self.r_s.append(rs)

	def compute_demand(self, dt: float):
		"""
		Computation of the demand of the ramp exiting this service station at time instant k

		Note: Now takes into account value prescribed from traffic controller
		"""

		if self.k < self.delta:
			d_s = self.e[self.k] / dt
		else:
			d_s = self.s_s[self.k - round(self.delta)] + (self.e[self.k] / dt)

		self.dem.append(min([d_s, self.r_s_c[self.k, self.id], self.r_s_max]))

	def compute_num_vehicles(self, dt: float):
		"""
		Computation of the number of vehicles at this service station at time instant k + 1

		Note: Doesn't include vehicles in service-station queue
		"""

		if len(self.s_s) < self.delta:
			self.l.append(self.l[self.k] + dt * self.s_s[self.k])
		else:
			self.l.append(self.l[self.k] + dt * (self.s_s[self.k] - self.s_s[int(self.k - self.delta)]))

	def compute_queue_length(self, dt: float):
		"""
		Computation of the number of vehicles queueing at this service station at time instant k + 1
		(due to the impossibility of merging back into the mainstream)
		"""

		d_queue = self.e[self.k] - (dt * self.r_s[-1])
		s_queue = 0

		if len(self.s_s) >= self.delta:
			s_queue += (dt * self.s_s[int(self.k - self.delta)])
			d_queue += s_queue

		self.s_e.append(s_queue)
		self.e.append(d_queue)

	def update_k(self, kappa: int):
		"""
		Each iteration starts with the update of the time instant
		"""

		self.k = kappa

	def reset(self):
		"""
		Reset station to zero initial condition
		"""

		self.k = 0
		self.s_s = []
		self.s_e = []
		self.r_s = []
		self.e = [0]
		self.l = [0]
		self.dem = []
