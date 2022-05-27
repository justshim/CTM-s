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
		## Utility method to print some information about the cell

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
		## Computation of the flow entering this cell from the previous one at time instant k
		## (the computation of phi varies according to the congestion state of the cell, hence we update it first)
 
		self.updateCongestionState(Dprec, TotalDs)

		if(self.congestionState == 0): #FREE FLOW
			self.phi = Dprec
			#print("Free flow, cell " + str(self.ID_cell))
			#print("Compute phi: Dprec " + str(Dprec))
		
		elif(self.congestionState == 1): #CONGESTED MAINSTREAM
			self.phi = self.SBig-TotalDs
			print("Congested 1, cell " + str(self.ID_cell))
			#print("TotalDs " + str(TotalDs))
			#print("SBig " + str(self.SBig))
		
		elif(self.congestionState == 2): #CONGESTED SERVICE
			self.phi = Dprec
			print("Congested 2, cell " + str(self.ID_cell))
		
		elif(self.congestionState == 3): #CONGESTED ALL
			self.phi = self.SBig  * self.p_ms
			print("Congested 3, cell " + str(self.ID_cell))
		
		#print("Phi: " + str(self.phi))

	def computePhiPlus(self, Rs_total):
		## Computation of the total flow entering this cell at time instant k

		self.phi_plus = self.phi + self.r + Rs_total
		#print("phi +: " + str(self.phi_plus)) 
		
	def computePhiMinus(self, Ss, NextPhi):
		## Computation of the total flow exiting this cell at time instant k

		self.phi_minus = NextPhi + self.s + Ss
		#print("phi -: " + str(self.phi_minus))

	def computeRho(self, TimeLength):
		## Computation of the traffic density of this cell at time instant k + 1

		self.rho.append(self.rho[self.k] + (TimeLength/self.length * (self.phi_plus - self.phi_minus)))
		#print("rho[k+1]:  " + str(self.rho[self.k+1]))


	def computeDBig(self, total_beta):
		## Computation of the demand of this cell at time instant k

		supp = (1 - self.beta - total_beta) * self.v * self.rho[self.k]	# this is a support variable used to simplify the syntax later

		if(supp > self.q_max):
			self.DBig = self.q_max
		
		else:
			self.DBig = supp
		
		#print("DBig:  " + str(self.DBig))

	def computeSBig(self):
		## Computation of the supply of this cell at time instant k

		supp = self.w * (self.rho_max-self.rho[self.k]) 	# this is a support variable used to simplify the syntax later
		
		if(supp > self.q_max):
			self.SBig = self.q_max
		
		else:
			self.SBig = supp

	def updateCongestionState(self, Dprec, TotalDs):
		## Computation of the congestion state of this cell at time instant k

		#print("TotalDs: "+ str(TotalDs))
		if(Dprec + TotalDs <= self.SBig): 
			self.congestionState=0 #FREE FLOW
		
		elif((Dprec > self.p_ms * self.SBig) and (TotalDs <= (1 - self.p_ms) * self.SBig)):
			self.congestionState=1 #CONGESTED MAINSTREAM
		
		elif((Dprec <= self.p_ms * self.SBig) and (TotalDs > (1 - self.p_ms) * self.SBig)):
			self.congestionState=2 #CONGESTED SERVICE
		
		elif((Dprec > self.p_ms * self.SBig) and (TotalDs > (1 - self.p_ms) * self.SBig)):
			self.congestionState=3 #CONGESTED ALL
		
		else:
			print("Congestion Error")


	def updateK(self, kappa):
		## Each iteration starts with the update of the time instant

		self.k=kappa