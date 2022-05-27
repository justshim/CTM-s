class OnRamp:
	"""Class modeling the on-ramps on an highway stretch in the CTM-s model"""
	def __init__(self, ID, d_r, r_r_max, j, p_r):
		self.ID_onramp = ID
		self.j = j
		self.p_r = p_r
		self.r_r_max = r_r_max
		self.r_r = 0
		self.l_r = [0]
		self.d_r = d_r
		self.d_r_big = 0
		self.k = 0

	def toString(self):
		## Utility method to print some information about the on-ramp

		print("On-Ramp ID: "+str(self.ID_onramp))
		print("To cell "+str(self.j))
		print("Number of vehicles: "+str(self.l_r))
		print()

	def computeDrBig(self, timeLength):
		## Computation of the demand of this ramp at time instant k

		supp = self.d_r + self.l_r[self.k]/timeLength

		if(supp > self.r_r_max):
			self.d_r_big = self.r_r_max

		else:
			self.d_r_big = supp

	def computeL(self, timeLength):
		## Computation of the number of vehicles on this ramp at time instant k + 1
		
		self.l_r.append(self.l_r[self.k] + timeLength * (self.d_r - self.r_r))

	def computeRr(self, t, rr):
		## Computation of the flow merging into the mainstream from this ramp at time instant k
		
		if t == 0 or t == 1:
			self.r_r = self.d_r_big

		elif t == 2 or t == 3:
			self.r_r = rr

	def updateK(self, kappa):
		## Each iteration starts with the update of the time instant

		self.k = kappa
		