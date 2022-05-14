import cell as c
import station as st

class Stretch:
	"""Controller class for the model, represents the system at large"""

	def __init__(self, timeLength, lastPhi, first_DBig):
		self.cells = []
		self.stations = []
		self.n_cells = 0
		self.n_stations = 0
		self.timeLength = timeLength
		self.TTT = 0
		self.delta_big = 0
		self.pi = 0
		self.lastPhi = lastPhi
		self.first_DBig = first_DBig

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
		
	def createCell(self, length, v, w, q, q_max, s, r, rho_max, beta, p):
		cell = c.Cell(self.n_cells, length, v, w, q, q_max, s, r, rho_max, beta, p) 
		self.cells.append(cell)
		self.n_cells = self.n_cells + 1

	def createStation(self, r_s_max, i, j, delta, beta_s, p):
		station = st.Station(self.n_stations, r_s_max, i, j, delta, beta_s, p) 
		self.stations.append(station)
		self.n_stations = self.n_stations + 1

	def setT(self, newT):
		self.timeLength=newT

	def update(self):
		##Station update
		#cong = 0
		#for s in self.stations:
		#	for c in self.cells:
		#			if(s.i==c.ID_cell):
		#				s.computeSs(c.phi_minus)
		#			if(s.j==c.ID_cell):
		#				cong=c.congestionState
		#
		#	s.computeRs(cong)
		#	s.computeL(self.timeLength)
		#	s.computeE(self.timeLength)
		#	s.computeDsBig(self.timeLength)


		#Cell update
		cong = 0
		beta=0
		total_beta=0
		total_Rs = 0 
		total_Ds=0
		next_phi = 0
		prev_DBig = 0

		for i in range (len(self.cells)):
			
			# special treatment for last cell
			if((i+1) < (len(self.cells))):
				next_phi = self.cells[i+1].phi
			else:
				next_phi = self.lastPhi 		# Static assignment from data for last cell

			# special treatment for first cell
			if(i != 0):
				prev_DBig = self.cells[i-1].DBig
			else:
				prev_DBig = self.first_DBig		# Static assignment from data for first cell

			# first total_Ds needs to be updated for the computation of the congestion state
			for s in self.stations:
				if(s.j==i):
					total_Ds=total_Ds+s.d_s_big
					#### ATTENZIONE: Ds DA RIFERIRE A USCITA DELLA STAZIONE?

			self.cells[i].computePhi(prev_DBig, total_Ds) # calls cell.updateCongestionState

			# if cell has stations entering or exiting, those stations are updated
			for s in self.stations:
				if(s.i==i):
					s.computeSs(self.cells[i].phi_minus)
					ss=s.Ss
					beta=s.beta_s
					total_beta = total_beta + s.beta_s
				else:
					ss=0

				if(s.j==i) or (s.i==i):
					if self.cells[i].congestionState == 0 or self.cells[i].congestionState == 1:
						s.computeRs(self.cells[i].congestionState)
					elif self.cells[i].congestionState == 2:
						self.iterativeProcedure(i)
					elif self.cells[i].congestionState == 3:
						#self.iterativeProcedure2(i)
						pass

					s.computeL(self.timeLength)
					s.computeE(self.timeLength)
					s.computeDsBig(self.timeLength)
				
				if(s.j==i):
					total_Rs = total_Rs + s.Rs			

			self.cells[i].computeRho(self.timeLength, ss, next_phi, total_Rs)
			self.cells[i].computeDBig(total_beta)
			self.cells[i].computeSBig()
					

		def iterativeProcedure(self, i):
			demands = [] # contains whole stations for convenience
			Rs_vector = []
			prev_D = self.cells[i-1].DBig
			supply = self.cells[i].SBig
			supply_cap = supply
			good = [0]
			sum_D_good = 0
			tol = 0.1 ### CHIEDERE QUALE TOLLERANZA USARE
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
						Rs_vector.append((d.ID_station, d.d_s_big))
						sum_D_good = sum_D_good + d.d_s_big
						supply_cap = supply_cap - d.d_s_big
				
			# compute sum of priorities for all involved stations
			for b in bad:
				sum_p = sum_p + b.p

			# compute remaining Rs
			for b in bad:
				Rs_vector.append((b.ID_station, (b.p/sum_p)*supply_cap))

			# update all Rs of all stations involved
			for k in range(len(Rs_vector)):
				for station in self.stations:
					if Rs_vector[k(0)] == station.ID_station:
						station.Rs = Rs_vector[k(1)]









					
			

		


