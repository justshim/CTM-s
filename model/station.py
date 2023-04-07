class Station:
	# TODO: Provide type hints for clarity
	# TODO: Clean up method function names
	# id_station: int 		# TODO: ...
	# r_s_max
	# r_s_c
	# i: int
	# j: int
	# s_s: List
	# r_s = 0
	# e = [0]
	e_max: int 				# Max Queue Length # TODO: ?
							# Max Queue Length can be determined by dividing length of ramp by average length of car...?
	# d_s_big = 0
	# delta = delta
	# beta_s = beta_s
	# l = [0]
	# p = p
	# k = 0

	"""Class modeling the service stations on a highway stretch in the CTM-s model"""
	def __init__(self, id_station, r_s_max, i, j, delta, beta_s, p):
		self.id_station = id_station
		self.r_s_max = r_s_max
		self.r_s_c = r_s_max
		self.i = i
		self.j = j
		self.delta = delta
		self.beta_s = beta_s
		self.e_max = 1000  # TODO: !!!
		self.p = p

		self.s_s = []
		self.r_s = []
		self.e = [0]
		self.l = [0]
		self.demand = 0
		self.k = 0

	def to_string(self):
		"""
		Utility method to print some information about the service station
		"""

		print("Station ID: " + str(self.id_station))
		print("From cell " + str(self.i))
		print("To cell " + str(self.j))
		print("Time delay: " + str(self.delta))
		print("Split ratio: " + str(self.beta_s))
		print("Number of vehicles: " + str(self.l))
		print()

	def compute_ss(self, next_phi):
		"""
		Computation of the flow leaving the mainstream to enter this service station at time instant k
		"""

		self.s_s.append((self.beta_s / (1 - self.beta_s)) * next_phi)
		
	def compute_rs(self, rs, t):
		"""
		Computation of the flow merging into the mainstream from this service station at time instant k
		"""

		if t == 0 or t == 1:
			self.r_s.append(self.demand)

		elif t == 2 or t == 3:
			self.r_s.append(rs)

	def compute_demand(self, dt):
		"""
		Computation of the demand of the ramp exiting this service station at time instant k
		Update: Now takes into account value prescribed from traffic controller
		In uncontrolled case,
		"""

		if self.k < self.delta:
			d_s = (0 + self.e[self.k]) / dt  # TODO: Check this???
		else:
			d_s = self.s_s[self.k - round(self.delta)] + self.e[self.k] / dt  # TODO: Check this??

		# TODO: Need to change this to r_s_c[k, id]
		self.demand = min([d_s, self.r_s_c, self.r_s_max])

	def compute_num_vehicles(self, dt: float):
		"""
		Computation of the number of vehicles at this service station at time instant k + 1
		"""

		## L and e together:
		# self.l.append(self.l[self.k] + time_length * self.s_s[self.k] - time_length * self.r_s)

		## L and e separated:
		if len(self.s_s) < self.delta:
			self.l.append(self.l[self.k] + dt * self.s_s[self.k] - 0)
		else:
			self.l.append(self.l[self.k] + dt * (self.s_s[self.k] - self.s_s[int(self.k - self.delta)]))

	def compute_queue_length(self, dt):
		"""
		Computation of the number of vehicles queueing at this service station at time instant k
		(due to the impossibility of merging back into the mainstream)
		"""

		d_queue = self.e[self.k] - (dt * self.r_s[-1])

		if len(self.s_s) >= self.delta:
			d_queue += (dt * self.s_s[int(self.k - self.delta)])

		self.e.append(d_queue)

	def update_k(self, kappa):
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
		self.r_s = []
		self.e = [0]
		self.demand = 0
		self.l = [0]
