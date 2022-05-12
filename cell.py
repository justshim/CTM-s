class Cell:
	"""Class describing the cells of the CTM model"""
	def __init__(self, ID, length, v, w, q, q_max, s, r, rho_max, beta, p):
		self.ID = ID
		self.length = length
		self.v = v
		self.w = w
		self.q = q
		self.p = p
		self.q_max = q_max
		self.s = s
		self.r = r
		self.rho_max = rho_max
		self.phi_minus = 0
		self.phi_plus = 0
		self.rho = 0
		self.beta = beta
		self.DBig = 0 
		self.SBig = 0
		self.congestionState = 0


	def toString(self):
		print("Cell ID: "+str(self.ID))
		print("Length: "+str(self.length))
		print("r: "+str(self.r))
		print("s: "+str(self.s))

	def computePhi(self, ID):
		return 0

	def computePhiPlus(self, Rs):
		#DA RIVEDERE PER L'1 TRASPOSTO
		self.phi_plus=self.computePhi(self.ID)+self.r+Rs
		return self.phi_plus

	def computePhiMinus(self, Ss):
		self.phi_minus=self.computePhi(int(self.ID)+1)+self.s+Ss
		return self.phi_minus

	def computeRho(self, TimeLength, Ss):
		phi_meno=self.computePhiMinus(Ss)
		phi_piu=self.computePhiPlus(Ss)
		self.rho=self.rho+(TimeLength/self.length*(phi_piu-phi_meno))
		

	def computeDBig(self, betaStation):
		a = (1-self.beta-betaStation)*self.v*self.rho
		if(a > self.q_max):
			self.DBig=self.q_max
		else:
			self.DBig=a

	def computeSBig(self):
		a = self.w*(self.rho_max-self.rho)
		if(a > self.q_max):
			self.SBig=self.q_max
		else:
			self.SBig=a

	def updateCongestionState(self, Dprec, TotalDs):
		if(Dprec+TotalDs<=self.SBig): 
			self.congestionState=0 #FREE FLOW
		#CONTROLLARE PMS O SOLO P
		elif((Dprec > self.p*self.SBig) and (TotalDs <= (1-self.p))):
			self.congestionState=1 #CONGESTED MAINSTREAM
		elif((Dprec <= self.p*self.SBig) and (TotalDs > (1-self.p))):
			self.congestionState=2 #CONGESTED SERVICE
		elif((Dprec > self.p*self.SBig) and (TotalDs > (1-self.p))):
			self.congestionState=3 #CONGESTED ALL
		else:
			print("Congestion Error")