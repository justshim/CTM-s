class Cell:
	# TODO: Provide type hints for clarity
	# TODO: Clean up method function names
	# TODO: Congestion state variables

	"""
	Class describing the cells of the CTM-s model
	"""
	def __init__(self, id_cell, length, v_free, w, q_max, rho_max, p_ms):
		self.id_cell = id_cell
		self.l = length
		self.v_free = v_free
		self.v = []
		self.w = w
		# self.q = [4000]  # Q is a vector for capacity drop modelling
		self.p_ms = p_ms
		self.q_max = q_max
		self.rho_max = rho_max

		self.rho = [0]
		self.phi = []
		self.phi_minus = 0
		self.phi_plus = 0
		self.demand = 0
		self.supply = 0
		self.cong_state = 0
		self.k = 0

	def to_string(self):
		"""
		Utility method to print some information about the cell
		"""

		print("Cell ID: " + str(self.id_cell))
		print("Length: " + str(self.l))
		print("v_free: " + str(self.v_free))
		print("w: " + str(self.w))
		print("rho_max: " + str(self.rho_max))
		print("q_max: " + str(self.q_max))
		print()

	# def computeQ(self):
	# 	# Q is a vector for capacity drop modelling
	# 	pass

	# def computeV(self):
	# 	#if(self.rho[self.k] == 0 or self.congestion_state == 0 or self.congestion_state == 2):
	# 	if(self.rho[self.k] == 0):
	#
	# 		self.v.append(self.v_free)
	# 	else:
	#
	# 		phi_avg = (self.phi_minus + self.phi_plus)/2
	# 		#print("rho: ")
	# 		#print(self.rho[self.k])
	# 		self.v.append(phi_avg/self.rho[self.k])

	def compute_phi(self, d_prev: float, total_ds: float):
		"""
		Computation of the flow entering this cell from the previous one at time instant k
		Varies according to the congestion state of the cell so we update it first
		"""

		self.update_cong_state(d_prev, total_ds)

		if self.cong_state == 0:  # FREE FLOW
			self.phi.append(d_prev)
		
		elif self.cong_state == 1:  # CONGESTED MAINSTREAM
			self.phi.append(self.supply - total_ds)

		elif self.cong_state == 2:  # CONGESTED SERVICE
			self.phi.append(d_prev)

		elif self.cong_state == 3:  # CONGESTED ALL
			self.phi.append(self.supply * self.p_ms)

	def compute_phi_plus(self, rs_total: float):
		"""
		Computation of the total flow entering this cell at time instant k
		"""

		self.phi_plus = self.phi[-1] + rs_total
		
	def compute_phi_minus(self, ss: float, next_phi: float):
		"""
		Computation of the total flow exiting this cell at time instant k
		"""

		self.phi_minus = next_phi + ss

	def compute_rho(self, dt: float):
		"""
		Computation of the traffic density of this cell at time instant k + 1
		"""

		self.rho.append(self.rho[self.k] + (dt/self.l) * (self.phi_plus - self.phi_minus))

	def compute_demand(self, total_beta: float):
		"""
		Computation of the demand of this cell at time instant k
		"""

		self.demand = min([(1 - total_beta) * self.v_free * self.rho[self.k], self.q_max])
		
	def compute_supply(self):
		"""
		Computation of the supply of this cell at time instant k
		"""

		self.supply = min([self.w * (self.rho_max - self.rho[self.k]), self.q_max])

	def update_cong_state(self, d_prev: float, total_ds: float):
		"""
		Computation of the congestion state of this cell at time instant k
		"""

		if d_prev + total_ds <= self.supply:
			self.cong_state = 0  # FREE FLOW
		
		elif (d_prev > self.p_ms * self.supply) and (total_ds <= (1 - self.p_ms) * self.supply):
			self.cong_state = 1  # CONGESTED MAINSTREAM

		elif (d_prev <= self.p_ms * self.supply) and (total_ds > (1 - self.p_ms) * self.supply):
			self.cong_state = 2  # CONGESTED SERVICE

		elif (d_prev > self.p_ms * self.supply) and (total_ds > (1 - self.p_ms) * self.supply):
			self.cong_state = 3  # CONGESTED ALL
		else:
			print("cell: " + str(self.id_cell) + " Congestion Error")

	def update_k(self, kappa):
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
		self.phi_minus = 0
		self.phi_plus = 0
		self.rho = [0]
		self.demand = 0
		self.supply = 0
		self.cong_state = 0
		self.k = 0
