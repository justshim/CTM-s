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
		pass

	def computePi(self):
		pass
		
	def createCell(self, length, v, w, q_max, s, r, rho_max, beta, p):
		cell = c.Cell(self.n_cells, length, v, w, q_max, s, r, rho_max, beta, p) 
		self.cells.append(cell)
		self.n_cells = self.n_cells + 1

	def createStation(self, r_s_max, i, j, delta, beta_s, p):
		station = st.Station(self.n_stations, r_s_max, i, j, delta, beta_s, p) 
		self.stations.append(station)
		self.n_stations = self.n_stations + 1

	def setT(self, newT):
		self.timeLength=newT

	def update(self, k):
		#print("Time instant: " + str(k))
		prev_DBig = 0
		totalDs = 0
		Ss_tot = 0

		## First of all update time instant for all cells
		for i in range (len(self.cells)):
			self.cells[i].updateK(k)

		## For stations, first update time instant and compute preliminary values
		for s in range (len(self.stations)):
			self.stations[s].updateK(k)
			self.stations[s].computeDsBig(self.timeLength)
			#self.stations[s].computeRs()
		
		## First stock of cell value updates, with special case for cell 0
		for i in range (len(self.cells)):
			#print("Cell: " + str(i))
			totalBeta = 0
			totalDs = 0
			
			## For each cell, check if any station stems from it, and sum all betas
			for s in range (len(self.stations)):
				if self.stations[s].i == i:
					totalBeta += self.stations[s].beta_s

			## For each cell, check if any station merges in it, and sum all Ds's
			for s in range (len(self.stations)):
				if self.stations[s].j == i:
					totalDs += self.stations[s].d_s_big

			self.cells[i].computeDBig(totalBeta)
			prev_DBig = 0
			
			# special treatment for first cell
			if(i != 0):

				prev_DBig = self.cells[i-1].DBig
				#print("prev_DBig: " + str(self.cells[i-1].DBig))
			else:
				prev_DBig = self.phi_zero[k]		## Static assignment from data for first cell
				#print("prev_DBig: " + str(prev_DBig))
			
			self.cells[i].computePhi(prev_DBig, totalDs)
			#print("Phi: " + str(self.cells[i].phi))

		## Second stock of cell value updates, with special case for cell 0
		for i in range (len(self.cells)):
			next_phi = 0
			totalRs = 0
			Ss_tot = 0
			#print("Cell: " + str(i))
			
			#special treatment for last cell
			if((i+1) < (len(self.cells))):
				next_phi = self.cells[i+1].phi
			
			else:
				next_phi = self.lastPhi	## Static assignment from data for last cell
			
			for s in range (len(self.stations)):
				
				if self.stations[s].j == i:
					#print("i: "+str(i)+"	self.stations[s].j: "+str(self.stations[s].j))
					if self.cells[i].congestionState == 0 or self.cells[i].congestionState == 1:
	 					self.stations[s].computeRs()
					
					elif self.cells[i].congestionState == 2:
	 					self.iterativeProcedure(i, 2)
					
					elif self.cells[i].congestionState == 3:
	 					self.iterativeProcedure(i, 3)
					
					totalRs += self.stations[s].Rs
				
					print("Rs: "+str(self.stations[s].Rs))
				
				if self.stations[s].i == i:
					self.stations[s].computeSs(next_phi)
					Ss_tot += self.stations[s].Ss
		
			#print("next_phi: " + str(next_phi))

			self.cells[i].computeSBig()
			self.cells[i].computePhiMinus(Ss_tot, next_phi)
			print("Total RS: "+str(totalRs))
			self.cells[i].computePhiPlus(totalRs)
			self.cells[i].computeRho(self.timeLength)

		for s in range (len(self.stations)):
			self.stations[s].computeL(self.timeLength)
			self.stations[s].computeE(self.timeLength)


	def iterativeProcedure(self, i, t):
		print("Itero")
		demands = [] ## contains whole stations for convenience
		Rs_vector = []
		prev_D = self.cells[i-1].DBig
		supply = self.cells[i].SBig
		supply_res = supply
		good = [0]
		sum_D_good = 0
		sum_p = 0

		for s in self.stations:
			if(s.j==i):
			 	demands.append(s)

		bad = demands
			
		# compute E_cal_overline and E_cal_underline
		while len(good) != 0:
			good.clear()
			for d in demands:
				if d.d_s_big <= (supply - prev_D - sum_D_good)/len(bad):
					bad.remove(d)
					good.append(d)
					Rs_vector.append([d.ID_station, d.d_s_big])
					sum_D_good = sum_D_good + d.d_s_big
					if t == 2:
						supply_res = supply_res - d.d_s_big
					elif t == 3:
						supply_res = (1 - self.cells[i].p_ms)*supply_res - d.d_s_big
						## VERIFICARE CHE SIA DAVVERO Pms
				
		# compute sum of priorities for all involved stations
		for b in bad:
			sum_p = sum_p + b.p

		# compute remaining Rs
		for b in bad:
			Rs_vector.append((b.ID_station, (b.p/sum_p)*supply_res))

		# update all Rs of all stations involved
		for k in range(len(Rs_vector)):
			for station in self.stations:
				if Rs_vector[[k][0]] == station.ID_station:
					station.Rs = Rs_vector[[k][1]]

		









					
			

		


