from model import cell as c
from model import station as st
from model import on_ramp as onr
from model import off_ramp as offr

class Stretch:
	"""Controller class for the model, represents the system at large"""

	def __init__(self, time_length, last_phi, phi_zero):
		self.cells = []
		self.stations = []
		self.on_ramps = []
		self.off_ramps = []
		self.n_cells = 0
		self.n_stations = 0
		self.n_on_ramps = 0
		self.n_off_ramps = 0
		self.time_length = time_length ### T[h]
		self.ttt = 0
		self.delta_big = []
		self.last_phi = last_phi
		self.phi_zero = phi_zero
		self.k = 0

	def toString(self):
		## Utility method to print some information about the highway stretch

		for i in self.cells:
			i.toString()
		for s in self.stations:
			s.toString()
		print("N_cells: "+str(self.n_cells))
		print("N_stations: "+str(self.n_stations))
		print("Time Length: "+str(self.time_length))
		print("TTT: "+str(self.ttt))
		print("delta_big: "+str(self.delta_big))
		print("pi: "+str(self.pi))
		print()

	def computeTTT(self):
		total_ell = 0
		for i in range(len(self.cells)):
			total_ell +=  self.cells[i].length/self.cells[i].v_free

		self.ttt = total_ell*3600
		#print("TTT: " + str(self.ttt))
		return self.ttt
		

	def computeDelta(self): #viene negativo a causa del calcolo della velocità 
		total_ell = 0
		## Computation of the additional TTT (total travel time) due to congestions on this stretch
		for i in range(len(self.cells)-1):
			total_ell += 60 * ((self.cells[i].length/self.cells[i].v[self.k]) - (self.cells[i].length/self.cells[i].v_free))
			
			if (total_ell<0):
				total_ell = 0

		self.delta_big.append(total_ell)
	
	def createCell(self, length, v_free, w, q_max, rho_max, p):
		## Method to create an instance of the object Cell, and add it to this stretch

		cell = c.Cell(self.n_cells, length, v_free, w, q_max, rho_max, p) 
		self.cells.append(cell)
		self.n_cells = self.n_cells + 1

	def createStation(self, r_s_max, i, j, delta, beta_s, p):
		## Method to create an instance of the object Station, and add it to this stretch

		station = st.Station(self.n_stations, r_s_max, i, j, delta, beta_s, p) 
		self.stations.append(station)
		self.n_stations = self.n_stations + 1

	def createOnRamp(self, d_r, r_r_max, j, p_r):
		## Method to create an instance of the object Station, and add it to this stretch

		on_ramp = onr.OnRamp(self.n_on_ramps, d_r, r_r_max, j, p_r) 
		self.on_ramps.append(on_ramp)
		self.n_on_ramps = self.n_on_ramps + 1

	def createOffRamp(self, i, beta_r):
		## Method to create an instance of the object Station, and add it to this stretch

		off_ramp = offr.OffRamp(self.n_off_ramps, i, beta_r) 
		self.off_ramps.append(off_ramp)
		self.n_off_ramps = self.n_off_ramps + 1

	def update(self, kappa):
		## Main method of the calss: at each time instant k updates all the parameters of the cells and service stations on this stretch
		# initialization of support variables
		self.k = kappa
		total_beta = 0
		total_ds = 0
		prev_d_big = 0
		next_phi = 0
		total_rs = 0
		total_ss = 0
		total_rs_station = 0
		total_ss_station = 0
		total_rr_ramp = 0
		total_sr_ramp = 0

		self.preliminary_updates()
		
		## First batch of cell value updates, with special case for cell 0
		for i in range (len(self.cells)):
			
			total_beta = self.computeTotalBeta(i)
			total_ds = self.computeTotalDs(i)

			self.cells[i].computeDBig(total_beta)
			
			prev_d_big = self.computeDPrec(i)
			
			self.cells[i].computePhi(prev_d_big, total_ds)

		
		## Second batch of cell value updates, with special case for last cell
		for i in range (len(self.cells)):
			next_phi = self.computeNextPhi(i)
			
			total_rs_station = self.computeRsStation(i)
			total_ss_station = self.computeSsStation(i, next_phi)
			
			total_sr_ramp = self.computeSrRamp(i, next_phi)
			total_rr_ramp = self.computeRrRamp(i)
			
			total_rs = total_rr_ramp + total_rs_station
			total_ss = total_sr_ramp + total_ss_station

			self.cells[i].computeSBig()
			self.cells[i].computePhiMinus(total_ss, next_phi)
			self.cells[i].computePhiPlus(total_rs)
			self.cells[i].computeRho(self.time_length)
			
		self.finalUpdates()
		

	def iterativeProcedure(self, i, t):
		## Method called during the update procedure and used to assign r_s to all stations merging into the same cell in case of congestions of type 2 and 3

		demands = []		# initialization of support variables
		prev_d = self.cells[i-1].d_big
		supply = self.cells[i].s_big
		good = [0]			# list to contain "good" demands, i.e. the ones that do not saturate the flow
		sum_d_good = 0
		sum_p = 0

		for s in self.stations:
			if s.j == i:
			 	demands.append(s)

		if t == 2:
			supply_res = supply - prev_d

		elif t == 3:
			supply_res = (1 - self.cells[i].p_ms) * supply

		bad = demands		# list to contain "bad" demands, i.e. the ones that do saturate the flow
			
		# Recursively compute "bad" (E_cal_overline in the paper) and "good" (E_cal_underline in the paper)
		while len(good) != 0:
			good.clear()
			if t == 2:
				for d in demands:
					if d.d_s_big <= (supply_res - sum_d_good)/len(bad):
						bad.remove(d)
						good.append(d)
						
						#update RS for "good"
						for station in self.stations: 
							if d.id_station == int(station.id_station):
								station.computeRs(d.d_s_big, t)

						sum_d_good = sum_d_good + d.d_s_big

						supply_res = supply_res - d.d_s_big
			elif t == 3:
				for d in demands:
					if d.d_s_big <= (((1 - self.cells[i].p_ms) * supply_res) - sum_d_good)/len(bad):
						bad.remove(d)
						good.append(d)
						
						#update RS for "good"
						for station in self.stations: 
							if d.id_station == int(station.id_station):
								station.computeRs(d.d_s_big, t)
						
						sum_D_good = sum_D_good + d.d_s_big

						supply_res = (1 - self.cells[i].p_ms) * supply_res - d.d_s_big  ## VERIFICARE FORMULA

		# Compute sum of priorities for all involved stations
		for b in bad:
			sum_p = sum_p + b.p

		# Update RS for "bad"
		for b in bad:	
			for station in self.stations: 
				if b.id_station == int(station.id_station):
					station.computeRs((b.p/sum_p) * supply_res, t)

	def preliminary_updates(self): 
		## First of all update time instant for all cells with current k
		for i in range (len(self.cells)):
			self.cells[i].updateK(self.k)

		## Same for stations, plus computation of some preliminary values
		for s in range (len(self.stations)):
			self.stations[s].updateK(self.k)
			self.stations[s].computeDsBig(self.time_length)

		## Same for on-ramps, plus computation of some preliminary values
		for r_on in range (len(self.on_ramps)):
			self.on_ramps[r_on].updateK(self.k)
			self.on_ramps[r_on].computeDrBig(self.time_length)

	def computeTotalBeta(self, i):
		total_beta = 0
		## For each cell, check if any station stems from it, and sum all betas (needed for the computation of D_i)
		for s in range (len(self.stations)):
			if self.stations[s].i == i:
				total_beta += self.stations[s].beta_s

		## For each cell, check if any station stems from it, and sum all betas (needed for the computation of D_i)
		for r_off in range (len(self.off_ramps)):
			if self.off_ramps[r_off].i == i:
				total_beta += self.off_ramps[r_off].beta_r
		
		return total_beta

	def computeTotalDs(self, i):
		total_ds = 0
		## For each cell, check if any station merges in it, and sum all Ds's (needed for the computation of phi_i)
		for s in range (len(self.stations)):
			if self.stations[s].j == i:
				total_ds += self.stations[s].d_s_big

		## For each cell, check if any on-ramp merges in it, and sum all Dr's (needed for the computation of phi_i)
		for r_on in range (len(self.on_ramps)):
			if self.on_ramps[r_on].j == i:
				total_ds += self.on_ramps[s].d_r_big

		return total_ds
	
	def computeRsStation(self, i):
		## For each cell, check if any stations merge into it, and compute their r_s; 
		total_rs = 0
		for s in range (len(self.stations)):	
			if self.stations[s].j == i:
				#print("i: "+str(i)+"	self.stations[s].j: "+str(self.stations[s].j))
				if self.cells[i].congestion_state == 0 or self.cells[i].congestion_state == 1:
	 				self.stations[s].computeRs(0, self.cells[i].congestion_state)
					
				elif self.cells[i].congestion_state == 2:
	 				self.iterativeProcedure(i, self.cells[i].congestion_state)
					
				elif self.cells[i].congestion_state == 3:
	 				self.iterativeProcedure(i, self.cells[i].congestion_state)
					
				total_rs += self.stations[s].r_s
				
		return total_rs

	def computeSsStation(self, i, next_phi):
		##check if any stations stem from it, and compute their s_s. These are then summed up for use, respectively, in the computation of Phi- and Phi+
		total_ss=0
		for s in range (len(self.stations)):
			if self.stations[s].i == i:
				self.stations[s].computeSs(next_phi)
				total_ss += self.stations[s].s_s[self.k]

		return total_ss

	def computeRrRamp(self, i):
		total_rs=0
		for r_on in range (len(self.on_ramps)):

			if self.on_ramps[r_on].j == i:
					
				if self.cells[i].congestion_state == 0 or self.cells[i].congestion_state == 1:
	 				self.on_ramps[r_on].computeRr(0, self.cells[i].congestion_state)
				
				elif self.cells[i].congestion_state == 2:
	 				rr = self.cells[i].s_big - self.cells[i-1].d_big
	 				self.on_ramps[r_on].computeRr(rr, self.cells[i].congestion_state)
				
				elif self.cells[i].congestion_state == 3:
	 				rr = self.cells[i].s_big * self.on_ramps[r_on].p_r
	 				self.on_ramps[r_on].computeRr(rr, self.cells[i].congestion_state)
				
				total_rs += self.on_ramps[r_on].r_r

		return total_rs

	def computeSrRamp(self, i, next_phi):
		total_ss=0
		for r_off in range (len(self.off_ramps)):

			if self.off_ramps[r_off].i == i:
				self.off_ramps[r_off].computeSr(next_phi)
				total_ss += self.off_ramps[r_off].s_r

		return total_ss

	def finalUpdates(self):
		## As a final step, all stations have their l and e updated
		for s in range (len(self.stations)):
			self.stations[s].computeE(self.time_length)
			self.stations[s].computeL(self.time_length)

		## And all ramps have their l updated
		for r_on in range (len(self.on_ramps)):
			self.on_ramps[r_on].computeL(self.time_length)

		for i in range(len(self.cells)):
			self.cells[i].computeV()
			

		self.computeDelta()

	def computeDPrec(self, i):
		prev_d_big = 0
		# First cell does not have a "previous" cell, hence phi_(i-1) is given as input
		if(i != 0):
			prev_d_big = self.cells[i-1].d_big
		
		else:
			prev_d_big = self.phi_zero[self.k]		
			
		return prev_d_big

	def computeNextPhi(self, i):
		next_phi=0
		# Last cell does not have a "next" cell, hence phi_(i+1) is given as input
		if((i+1) < (len(self.cells))):
			next_phi = self.cells[i+1].phi
		
		else:
			#next_phi = self.last_phi[self.k]
			next_phi = self.cells[i-1].d_big

		return next_phi
