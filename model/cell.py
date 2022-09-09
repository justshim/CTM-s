class Cell:
	"""Class describing the cells of the CTM model"""
	def __init__(self, id_cell, length, v_free, w, q_max, rho_max, p_ms):
		self.id_cell = id_cell
		self.length = length
		self.v_free = v_free
		self.v = []
		self.w = w
		self.q = [4000] # Q is a vector for capacity drop modelling
		self.p_ms = p_ms
		self.q_max = q_max
		self.rho_max = rho_max
		self.phi = 0
		self.phi_minus = 0
		self.phi_plus = 0
		self.rho = [0]
		self.d_big = 0 
		self.s_big = 0
		self.congestion_state = 0
		self.k = 0

	def toString(self):
		## Utility method to print some information about the cell

		print("Cell ID: "+str(self.id_cell))
		print("Length: "+str(self.length))
		print("v_free: "+str(self.v_free))
		print("w: "+str(self.w))
		print("rho_max: "+str(self.rho_max))
		print("q_max: "+str(self.q_max))
		print()

	def computeQ(self):
		#Q is a vector for capacity drop modelling
		pass

	def computeV(self):
		#if(self.rho[self.k] == 0 or self.congestion_state == 0 or self.congestion_state == 2):
		if(self.rho[self.k] == 0):
		
			self.v.append(self.v_free)
		else:

			phi_avg = (self.phi_minus + self.phi_plus)/2
			#print("phi_avg: ")
			#print(phi_avg)
			self.v.append(phi_avg/self.rho[self.k])

		

	def computePhi(self, d_prec, total_ds):
		## Computation of the flow entering this cell from the previous one at time instant k
		## (the computation of phi varies according to the congestion state of the cell, hence we update it first)
 
		self.updateCongestionState(d_prec, total_ds)

		if(self.congestion_state == 0): #FREE FLOW
			self.phi = d_prec
		
		elif(self.congestion_state == 1): #CONGESTED MAINSTREAM
			self.phi = self.s_big - total_ds
			#print("Congested 1, cell " + str(self.id_cell))

		elif(self.congestion_state == 2): #CONGESTED SERVICE
			self.phi = d_prec
			#print("Congested 2, cell " + str(self.id_cell))
		
		elif(self.congestion_state == 3): #CONGESTED ALL
			self.phi = self.s_big  * self.p_ms
			#print("Congested 3, cell " + str(self.id_cell))

		#print("congestion_state: ")
		#print(self.congestion_state)

	def computePhiPlus(self, rs_total):
		## Computation of the total flow entering this cell at time instant k

		self.phi_plus = self.phi + rs_total 
		
	def computePhiMinus(self, ss, next_phi):
		## Computation of the total flow exiting this cell at time instant k

		self.phi_minus = next_phi + ss

	def computeRho(self, time_length):
		## Computation of the traffic density of this cell at time instant k + 1

		self.rho.append(self.rho[self.k] + (time_length/self.length * (self.phi_plus - self.phi_minus)))

	def computeDBig(self, total_beta):
		## Computation of the demand of this cell at time instant k

		supp = (1 - total_beta) * self.v_free * self.rho[self.k]	# this is a support variable used to simplify the syntax later
		
		if(supp > self.q_max):
			self.d_big = self.q_max
		
		else:
			self.d_big = supp
		
	def computeSBig(self):
		## Computation of the supply of this cell at time instant k

		supp = self.w * (self.rho_max-self.rho[self.k]) 	# this is a support variable used to simplify the syntax later
		#if(self.k>2500 and self.k<3000):
				#print("Cell " + str(self.id_cell) + " k: " +str(self.k) + "\nsupp:" + str(supp) + " qmax: "+str(self.q_max))
		
		if(supp > self.q_max):	
			self.s_big = self.q_max
		else:
			self.s_big = supp

	def updateCongestionState(self, d_prec, total_ds):
		## Computation of the congestion state of this cell at time instant k

		if(d_prec + total_ds <= self.s_big): 
			self.congestion_state=0 #FREE FLOW
		
		elif((d_prec > self.p_ms * self.s_big) and (total_ds <= (1 - self.p_ms) * self.s_big)):
			self.congestion_state=1 #CONGESTED MAINSTREAM
			#print("Cell " + str(self.id_cell) + " in congestion")
		
		elif((d_prec <= self.p_ms * self.s_big) and (total_ds > (1 - self.p_ms) * self.s_big)):
			self.congestion_state=2 #CONGESTED SERVICE
			#print("Congestion ")
		
		elif((d_prec > self.p_ms * self.s_big) and (total_ds > (1 - self.p_ms) * self.s_big)):
			self.congestion_state=3 #CONGESTED ALL
			#print("Congestion ")
		
		else:
			print("cell: " + str(self.id_cell) + " Congestion Error")
			#print("D_prec: " + str(d_prec) + " self.p_ms * self.s_big: " + str(self.p_ms * self.s_big))
			#print("total_ds: " + str(total_ds) + " (1 - self.p_ms) * self.s_big): " + str((1 - self.p_ms) * self.s_big))


	def updateK(self, kappa):
		## Each iteration starts with the update of the time instant

		self.k=kappa