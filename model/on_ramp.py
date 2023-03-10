class OnRamp:
	# TODO: Provide type hints for clarity
	# TODO: Clean up method function names

	"""
	Class modeling the on-ramps on a highway stretch in the CTM-s model
	"""
	def __init__(self, id_onramp, d_r, r_r_max, j, p_r):
		self.id_onramp = id_onramp
		self.j = j
		self.p_r = p_r
		self.r_r_max = r_r_max
		self.r_r = 0
		self.l_r = [0]
		self.d_r = d_r
		self.d_r_big = 0
		self.k = 0

	def toString(self):
		"""
		Utility method to print some information about the on-ramp
		"""

		print("On-Ramp ID: "+str(self.id_onramp))
		print("To cell "+str(self.j))
		print("Number of vehicles: "+str(self.l_r))
		print()

	def computeDrBig(self, time_length):
		"""
		Computation of the demand of this ramp at time instant k
		"""

		supp = self.d_r + self.l_r[self.k]/time_length

		if supp > self.r_r_max:
			self.d_r_big = self.r_r_max

		else:
			self.d_r_big = supp

	def computeL(self, time_length):
		"""
		Computation of the number of vehicles on this ramp at time instant k + 1
		"""
		
		self.l_r.append(self.l_r[self.k] + time_length * (self.d_r - self.r_r))

	def computeRr(self, t, rr):
		"""
		Computation of the flow merging into the mainstream from this ramp at time instant k
		"""
		
		if t == 0 or t == 1:
			self.r_r = self.d_r_big

		elif t == 2 or t == 3:
			self.r_r = rr

	def updateK(self, kappa):
		"""
		Each iteration starts with the update of the time instant
		"""

		self.k = kappa
		