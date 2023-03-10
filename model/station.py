class Station:
	# TODO: Provide type hints for clarity
	# TODO: Clean up method function names

	"""Class modeling the service stations on a highway stretch in the CTM-s model"""
	def __init__(self, id_station, r_s_max, i, j, delta, beta_s, p):
		self.id_station = id_station
		self.r_s_max = r_s_max
		self.i = i
		self.j = j
		self.s_s = []
		self.r_s = 0
		self.e = [0]
		self.d_s_big = 0
		self.delta = delta
		self.beta_s = beta_s
		self.l = [0]
		self.p = p
		self.k = 0

	def toString(self):
		"""
		Utility method to print some information about the service station
		"""

		print("Station ID: "+str(self.id_station))
		print("From cell "+str(self.i))
		print("To cell "+str(self.j))
		print("Time delay: "+str(self.delta))
		print("Split ratio: "+str(self.beta_s))
		print("Number of vehicles: "+str(self.l))
		print()

	def computeSs(self, next_phi):
		"""
		Computation of the flow leaving the mainstream to enter this service station at time instant k
		"""

		self.s_s.append((self.beta_s / (1 - self.beta_s)) * next_phi)
		
	def computeRs(self, rs, t):
		"""
		Computation of the flow merging into the mainstream from this service station at time instant k
		"""

		if t == 0 or t == 1:
			self.r_s = self.d_s_big

		elif t == 2 or t == 3:
			self.r_s = rs

	def computeE(self, time_length):
		"""
		Computation of the number of vehicles queueing at this service station at time instant k
		(due to the impossibility of merging back into the mainstream)
		"""

		if len(self.s_s) < self.delta:
			self.e.append(self.e[self.k] + 0 - (time_length * self.r_s))

		else:
			self.e.append(self.e[self.k] + (time_length * self.s_s[self.k - self.delta]) - (time_length * self.r_s))

	def computeDsBig(self, time_length):
		"""
		Computation of the demand of the ramp exiting this service station at time instant k
		"""

		supp = 0
		
		if len(self.s_s) < self.delta:  # for the first delta time instants we skip the computation of s_s(k-delta), as it would send the index out of bounds, and s_s is zero in this period anyways
			supp = (0 + self.e[self.k]) / time_length
		else:
			supp = (self.s_s[self.k - self.delta] + self.e[self.k] / time_length)
			
		if supp > self.r_s_max:
			self.d_s_big = self.r_s_max
		else:
			self.d_s_big = supp

	def computeL(self, time_length):
		"""
		Computation of the number of vehicles at this service station at time instant k + 1
		"""
		
		## L and e together: 
		# self.l.append(self.l[self.k] + time_length * self.s_s[self.k] - time_length * self.r_s)
		
		## L and e separated: 
		if len(self.s_s) < self.delta: 
			self.l.append(self.l[self.k] + time_length * self.s_s[self.k] - 0)
		else:
			self.l.append(self.l[self.k] + time_length * (self.s_s[self.k] - self.s_s[self.k - self.delta]))

	def updateK(self, kappa):
		"""
		Each iteration starts with the update of the time instant
		"""
		self.k = kappa
		

