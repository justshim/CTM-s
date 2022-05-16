class Cell:
	"""Class describing the cells of the CTM model"""
	def __init__(self, ID, length, v, w, q, q_max, s, r, rho_max, beta, p_ms):
		self.ID_cell = ID
		self.length = length
		self.v = v
		self.w = w
		self.q = q
		self.p_ms = p_ms
		self.q_max = q_max
		self.s = s
		self.r = r
		self.rho_max = rho_max
		self.phi = 0
		self.phi_minus = 0
		self.phi_plus = 0
		self.rho = 0
		self.beta = beta
		self.DBig = 0 
		self.SBig = 0
		self.congestionState = 0


	def toString(self):
		print("Cell ID: "+str(self.ID_cell))
		print("Length: "+str(self.length))
		print("r: "+str(self.r))
		print("s: "+str(self.s))
		print("v: "+str(self.v))
		print("w: "+str(self.w))
		print("rho_max: "+str(self.rho_max))
		print("q_max: "+str(self.q_max))
		print()

	def computePhi(self, Dprec, TotalDs):
		self.updateCongestionState(Dprec, TotalDs)

		if(self.congestionState == 0): #FREE FLOW
			self.phi=Dprec
			print("Free flow")
		elif(self.congestionState == 1): #CONGESTED MAINSTREAM
			self.phi=self.SBig-TotalDs
			print("Congested 1")
		elif(self.congestionState == 2): #CONGESTED SERVICE
			self.phi=Dprec
			print("Congested 2")
		elif(self.congestionState == 3): #CONGESTED ALL
			self.phi=self.SBig*self.p_ms
			print("Congested 3")

		print("Phi " + str(self.phi))
 
	def computePhiPlus(self, Rs_total):
		self.phi_plus = self.phi + self.r + Rs_total
		print("phi + " + str(self.phi_plus)) 
		return self.phi_plus

	def computePhiMinus(self, Ss, NextPhi):
		self.phi_minus = NextPhi + self.s + Ss
		print("phi - " + str(self.phi_minus))
		return self.phi_minus

	def computeRho(self, TimeLength, Ss, NextPhi, Rs_total):
		phi_meno=self.computePhiMinus(Ss, NextPhi)
		phi_piu=self.computePhiPlus(Rs_total)
		self.rho=self.rho+(TimeLength/self.length*(phi_piu-phi_meno))


	def computeDBig(self, total_beta):
		a = (1 - self.beta - total_beta) * self.v * self.rho
		
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
		
		elif((Dprec > self.p_ms*self.SBig) and (TotalDs <= (1-self.p_ms))):
			self.congestionState=1 #CONGESTED MAINSTREAM
		
		elif((Dprec <= self.p_ms*self.SBig) and (TotalDs > (1-self.p_ms))):
			self.congestionState=2 #CONGESTED SERVICE
		
		elif((Dprec > self.p_ms*self.SBig) and (TotalDs > (1-self.p_ms))):
			self.congestionState=3 #CONGESTED ALL
		
		else:
			print("update")
			print(Dprec)
			#print(TotalDs)
			#print(self.SBig)
			#print("Congestion Error")