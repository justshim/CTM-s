class Cell(object):
	"""Class describing the cells of the CTM model"""
	def __init__(self, ID, length, v, w, q, q_max, phi, phi_plus, phi_minus, s, s_big, d_big, r, rho, rho_max, beta, p, is_congested):
		self.ID = ID
		self.length = length
		self.v = v
		self.w = w
		self.q = q
		self.q_max = q_max
		self.phi = phi
		self.phi_plus = phi_plus
		self.phi_minus = phi_minus
		self.s = s
		self.s_big = s_big
		self.d_big = d_big
		self.r = r
		self.rho = rho
		self.rho_max = rho_max
		self.beta = beta
		self.p = p
		self.is_congested = is_congested


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