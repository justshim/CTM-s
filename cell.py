class Cell:
	"""Class describing the cells of the CTM model"""
	def __init__(self, ID, length, v, w, q, q_max, s, r, rho_max, beta, p):
		self.ID = ID
		self.length = length
		self.v = v
		self.w = w
		self.q = q
		self.q_max = q_max
		self.s = s
		self.r = r
		self.rho_max = rho_max


	def toString(self):
		print("ID: "+str(self.ID))
		print("Length: "+str(self.length))

	def computePhi(self):
		pass

	def computePhiPlus(self):
		pass

	def computePhiMinus(self):
		pass

	def computeRho(self):
		pass

	def computeDBig(self):
		pass

	def computeSBig(self):
		pass

	def checkCongestion(self):
		pass