import cell as c
import station as st

class Stretch:
	"""Controller class for the model, represents the system at large"""

	def __init__(self, timeLength, lastPhi, phi_zero):
		self.cells = []
		self.stations = []
		self.n_cells = 0
		self.n_stations = 0
		self.timeLength = timeLength ### T[h]
		self.TTT = 0
		self.delta_big = 0
		self.pi = 0
		self.lastPhi = lastPhi
		self.phi_zero = phi_zero

	def toString(self):
		## Utility method to print some information about the highway stretch

		for c in self.cells:
			c.toString()
		for s in self.stations:
			s.toString()
		print("N_cells: "+str(self.n_cells))
		print("N_stations: "+str(self.n_stations))
		print("timeLength: "+str(self.timeLength))
		print("TTT: "+str(self.TTT))
		print("delta_big: "+str(self.delta_big))
		print("pi: "+str(self.pi))
		print()


	def computeDelta(self):
		## Computation of the additional TTT (total travel time) due to congestions on this stretch

		pass

	def computePi(self):
		## Computation of percentage of peak congestion reduction on this stretch

		pass
		
	def createCell(self, length, v, w, q_max, s, r, rho_max, beta, p):
		## Method to create an instance of the object Cell, and add it to this stretch

		cell = c.Cell(self.n_cells, length, v, w, q_max, s, r, rho_max, beta, p) 
		self.cells.append(cell)
		self.n_cells = self.n_cells + 1

	def createStation(self, r_s_max, i, j, delta, beta_s, p):
		## Method to create an instance of the object Station, and add it to this stretch

		station = st.Station(self.n_stations, r_s_max, i, j, delta, beta_s, p) 
		self.stations.append(station)
		self.n_stations = self.n_stations + 1

	def setT(self, newT):
		## Method to externally set the time instant length, mainly for debugging

		self.timeLength=newT

	def update(self, k):
		## Main method of the calss: at each time instant k updates all the parameters of the cells and service stations on this stretch

		#print("Time instant: " + str(k))
		prev_DBig = 0		# initialization of support variables used later
		totalDs = 0
		Ss_tot = 0

		## First of all update time instant for all cells with current k
		for i in range (len(self.cells)):
			self.cells[i].updateK(k)

		## Samefor stations, plus computation of some preliminary values
		for s in range (len(self.stations)):
			self.stations[s].updateK(k)
			self.stations[s].computeDsBig(self.timeLength)
			#self.stations[s].computeRs()
		
		## First batch of cell value updates, with special case for cell 0
		for i in range (len(self.cells)):
			#print("Cell: " + str(i))
			totalBeta = 0		# initialization of support variables
			totalDs = 0
			prev_DBig = 0

			## For each cell, check if any station stems from it, and sum all betas (needed for the computation of D_i)
			for s in range (len(self.stations)):
				if self.stations[s].i == i:
					totalBeta += self.stations[s].beta_s

			## For each cell, check if any station merges in it, and sum all Ds's (needed for the computation of phi_i)
			for s in range (len(self.stations)):
				if self.stations[s].j == i:
					totalDs += self.stations[s].d_s_big

			self.cells[i].computeDBig(totalBeta)
			
			# First cell does not have a "previous" cell, hence phi_(i-1) is given as input
			if(i != 0):
				prev_DBig = self.cells[i-1].DBig
				#print("prev_DBig: " + str(self.cells[i-1].DBig))
			
			else:
				prev_DBig = self.phi_zero[k]		
				#print("prev_DBig: " + str(prev_DBig))
			
			self.cells[i].computePhi(prev_DBig, totalDs)
			#print("Phi: " + str(self.cells[i].phi))

		## Second batch of cell value updates, with special case for last cell
		for i in range (len(self.cells)):
			next_phi = 0		# initialization of support variables
			totalRs = 0
			Ss_tot = 0
			#print("Cell: " + str(i))
			
			# Last cell does not have a "next" cell, hence phi_(i+1) is given as input
			if((i+1) < (len(self.cells))):
				next_phi = self.cells[i+1].phi
			
			else:
				next_phi = self.lastPhi
			
			## For each cell, check if any stations merge into it, and compute their r_s; 
			## then check if any stations stem from it, and compute their s_s. These are then summed up for use, respectively, in the computation of Phi- and Phi+
			for s in range (len(self.stations)):
				
				if self.stations[s].j == i:
					#print("i: "+str(i)+"	self.stations[s].j: "+str(self.stations[s].j))
					if self.cells[i].congestionState == 0 or self.cells[i].congestionState == 1:
	 					self.stations[s].computeRs()
					
					elif self.cells[i].congestionState == 2:
	 					self.iterativeProcedure(i, 2, k)
					
					elif self.cells[i].congestionState == 3:
	 					self.iterativeProcedure(i, 3, k)
					
					totalRs += self.stations[s].Rs[k]
				
					print("Rs: "+str(self.stations[s].Rs[k]))
				
				if self.stations[s].i == i:
					self.stations[s].computeSs(next_phi)
					Ss_tot += self.stations[s].Ss[k]
					
		
			#print("next_phi: " + str(next_phi))

			self.cells[i].computeSBig()
			self.cells[i].computePhiMinus(Ss_tot, next_phi)
			#print("Total RS: "+str(totalRs))
			self.cells[i].computePhiPlus(totalRs)
			self.cells[i].computeRho(self.timeLength)

		## As a final step, all stations have their l and e updated
		for s in range (len(self.stations)):
			self.stations[s].computeE(self.timeLength)
			self.stations[s].computeL(self.timeLength)


	def iterativeProcedure(self, i, t, k):
		## Method called during the update procedure and used to assign r_s to all stations merging into the same cell in case of congestions of type 2 and 3

		#print("Iterative procedure in process")

		demands = []		# initialization of support variables
		#Rs_vector = []
		prev_D = self.cells[i-1].DBig
		supply = self.cells[i].SBig
		good = [0]			# list to contain "good" demands, i.e. the ones that do not saturate the flow
		sum_D_good = 0
		sum_p = 0

		for s in self.stations:
			if(s.j==i):
			 	demands.append(s)

		if t == 2:
			supply_res = supply - prev_D

		elif t == 3:
			supply_res = (1 - self.cells[i].p_ms)*supply

		bad = demands		# list to contain "bad" demands, i.e. the ones that do saturate the flow
			
		# Recursively compute "bad" (E_cal_overline in the paper) and "good" (E_cal_underline in the paper)
		while len(good) != 0:
			good.clear()
			if t == 2:
				for d in demands:
					if d.d_s_big <= (supply_res - sum_D_good)/len(bad):
						bad.remove(d)
						good.append(d)
						
						#update RS for "good"
						for station in self.stations: 
							if d.ID_station == int(station.ID_station):
								station.Rs.append (d.d_s_big)

						sum_D_good = sum_D_good + d.d_s_big

						supply_res = supply_res - d.d_s_big
			elif t == 3:
				for d in demands:
					if d.d_s_big <= (((1 - self.cells[i].p_ms) * supply_res) - sum_D_good)/len(bad):
						bad.remove(d)
						good.append(d)
						
						#update RS for "good"
						for station in self.stations: 
							if d.ID_station == int(station.ID_station):
								station.Rs.append (d.d_s_big)
						
						sum_D_good = sum_D_good + d.d_s_big

						supply_res = (1 - self.cells[i].p_ms) * supply_res - d.d_s_big  ## VERIFICARE FORMULA

		# Compute sum of priorities for all involved stations
		for b in bad:
			sum_p = sum_p + b.p

		#update RS for "bad"
		for b in bad:	
			for station in self.stations: 
				if b.ID_station == int(station.ID_station):
					station.Rs.append((b.p/sum_p)*supply_res)
					#print("Stazione n: "+str(station.ID_station) + " RS:" + str(station.Rs[k]))
		


		# # Compute remaining Rs
		# for b in bad:
		# 	Rs_vector.append((b.ID_station, (b.p/sum_p)*supply_res))

		# # Update all Rs of all stations involved
		# for j in range(len(Rs_vector)):
		# 	for station in self.stations:
		# 		#print (station.ID_station)
		# 		if Rs_vector[j][0] == int(station.ID_station):
		# 			station.Rs = Rs_vector[j][1] 


					

		









					
			

		


