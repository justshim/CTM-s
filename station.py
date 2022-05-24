import supervisor as su

class Station:
	"""Class modeling the service stations on an highway stretch in the CTM-s model"""
	def __init__(self, ID, r_s_max, i, j, delta, beta_s, p):
		self.ID_station = ID
		self.r_s_max = r_s_max
		self.i = i
		self.j = j
		self.Ss = []
		self.Rs = []
		self.E = [0]
		self.d_s_big = 0
		self.delta = delta
		self.beta_s = beta_s
		self.l = [0]
		self.p = p
		self.k = 0

	def toString(self):
		## Utility method to print some information about the service station

		print("Station ID: "+str(self.ID_station))
		print("From cell "+str(self.i))
		print("To cell "+str(self.j))
		print("Time delay: "+str(self.delta))
		print("Split ratio: "+str(self.beta_s))
		print("Vehicles number: "+str(self.l))
		print()

	def computeSs(self, nextPhi):
		## Computation of the flow leaving the mainstream to enter this service station at time instant k

		self.Ss.append((self.beta_s / (1 - self.beta_s)) * nextPhi)
		print("Ss: "+str(self.Ss[self.k]))

	def computeRs(self):
		## Computation of the flow merging into the mainstream from this service station at time instant k

		self.Rs.append(self.d_s_big)
		#print("d_s_big " + str(self.d_s_big))

	def computeE(self, TimeLength):
		## Computation of the number of vehicles queueing at this service station at time instant k (due to the impossibility of merging back into the mainstream)

		if len(self.Ss) < self.delta:
			self.E.append(self.E[self.k] + 0 - (TimeLength * self.Rs[self.k]))
			print("IF")
			#self.E.append(0)
		else:
			self.E.append(self.E[self.k] + (TimeLength * self.Ss[self.k-self.delta]) - (TimeLength * self.Rs[self.k]))
			print("ELSE // Ss= " + str(self.Ss[self.k-self.delta]))

	def computeDsBig(self, TimeLength):
		## Computation of the demand of the ramp exiting this service station at time instant k

		app = 0		# this is a support variable used to simplify the syntax later
		
		if len(self.Ss) < self.delta:  # for the first delta time instants we skip the computation of s_s(k-delta), as it would send the index out of bounds, and s_s is zero in this period anyways
			app = (0 + self.E[self.k]) / TimeLength
			print("self.l[self.k]" + str(self.l[self.k]))
			print("self.E[self.k]" + str(self.E[self.k]))

		else:
			print("self.l[self.k]" + str(self.l[self.k]))
			print("self.E[self.k]" + str(self.E[self.k]))
			#print("self.Ss[self.k-self.delta]" + str(self.Ss[self.k-self.delta]))
			app = (self.Ss[self.k-self.delta] + self.E[self.k]/TimeLength)
			
		#print("app: "+str(app))
		
		if(app > self.r_s_max):
			self.d_s_big = self.r_s_max
		else:
			self.d_s_big = app

	def computeL(self, TimeLength):
		## Computation of the number of vehicles at this service station at time instant k + 1

		self.l.append(self.l[self.k] + (self.Ss[self.k]-self.Rs[self.k]))
		#print("l: " +str(self.l[self.k]))

	def updateK(self, kappa):
		## Each iteration starts with the update of the time instant
		
		self.k=kappa
		

