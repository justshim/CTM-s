class Cell:
	"""Class describing the cells of the CTM model"""
	def __init__(self, ID, length, v, w, q_max, s, r, rho_max, beta, p_ms):
		self.ID_cell = ID
		self.length = length
		self.v = v
		self.w = w
		self.q = [1500] # Q is a vector for capacity drop modelling
		self.p_ms = p_ms
		self.q_max = q_max
		self.s = s
		self.r = r
		self.rho_max = rho_max
		self.phi = 0
		self.phi_minus = 0
		self.phi_plus = 0
		self.rho = [0]
		self.beta = beta
		self.DBig = 0 
		self.SBig = 0
		self.congestionState = 0
		self.k = 0

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

	def computeQ(self):
		#Q is a vector for capacity drop modelling
		pass

	def computePhi(self, Dprec, TotalDs):
		self.updateCongestionState(Dprec, TotalDs)

		if(self.congestionState == 0): #FREE FLOW
			self.phi=Dprec
			print("Free flow")
			#print("Compute phi: Dprec " + str(Dprec))
		
		elif(self.congestionState == 1): #CONGESTED MAINSTREAM
			self.phi=self.SBig-TotalDs
			print("Congested 1")
			#print("TotalDs " + str(TotalDs))
			#print("SBig " + str(self.SBig))
		
		elif(self.congestionState == 2): #CONGESTED SERVICE
			self.phi=Dprec
			print("Congested 2")
		
		elif(self.congestionState == 3): #CONGESTED ALL
			self.phi=self.SBig*self.p_ms
			print("Congested 3")
		
		#print("Phi: " + str(self.phi))

	def computePhiPlus(self, Rs_total):
		self.phi_plus = self.phi + self.r + Rs_total
		#print("phi +: " + str(self.phi_plus)) 
		
	def computePhiMinus(self, Ss, NextPhi):
		self.phi_minus = NextPhi + self.s + Ss
		#print("phi -: " + str(self.phi_minus))

	def computeRho(self, TimeLength):
		self.rho.append(self.rho[self.k]+(TimeLength/self.length*(self.phi_plus-self.phi_minus)))
		#print("rho[k+1]:  " + str(self.rho[self.k+1]))


	def computeDBig(self, total_beta):
		a = (1 - self.beta - total_beta) * self.v * self.rho[self.k]

		if(a > self.q_max):
			self.DBig=self.q_max
		
		else:
			self.DBig=a
		
		#print("DBig:  " + str(self.DBig))

	def computeSBig(self):
		a = self.w*(self.rho_max-self.rho[self.k])
		
		if(a > self.q_max):
			self.SBig=self.q_max
		
		else:
			self.SBig=a

	def updateCongestionState(self, Dprec, TotalDs):
		#print("TotalDs: "+ str(TotalDs))
		if(Dprec+TotalDs<=self.SBig): 
			self.congestionState=0 #FREE FLOW
		
		elif((Dprec > self.p_ms*self.SBig) and (TotalDs <= (1-self.p_ms))):
			self.congestionState=1 #CONGESTED MAINSTREAM
		
		elif((Dprec <= self.p_ms*self.SBig) and (TotalDs > (1-self.p_ms))):
			self.congestionState=2 #CONGESTED SERVICE
		
		elif((Dprec > self.p_ms*self.SBig) and (TotalDs > (1-self.p_ms))):
			self.congestionState=3 #CONGESTED ALL
		
		else:
			print("Congestion Error")


	def updateK(self, kappa):
		self.k=kappa
		

	# def setPhiPlus(self, phi_uno):
	# 	self.phi_plus=phi_uno
	# 	#print("phi +: " + str(self.phi_plus)) 