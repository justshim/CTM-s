import supervisor as su

class Station:
	"""Class modeling the service stations on an highway stretch in the CTM-s model"""
	def __init__(self, ID, r_s_max, i, j, delta, beta_s, p):
		self.ID_station = ID
		self.r_s_max = r_s_max
		self.i = i
		self.j = j
		self.Ss = 0
		self.Rs = 0
		self.E = [0]
		self.d_s_big = 0
		self.delta = delta
		self.beta_s = beta_s
		self.l = [0]
		self.p = p
		self.k = 0

	def toString(self):
		print("Station ID: "+str(self.ID_station))
		print("From cell "+str(self.i))
		print("To cell "+str(self.j))
		print("Time delay: "+str(self.delta))
		print("Split ratio: "+str(self.beta_s))
		print("Vehicles number: "+str(self.l))
		print()

	def computeSs(self, nextPhi):
		self.Ss = (self.beta_s / (1 - self.beta_s)) * nextPhi ### ATTENZIONE: SOMMA DEI BETA?

	def computeRs(self):
			self.Rs=self.d_s_big
			#print("d_s_big " + str(self.d_s_big))

	def computeE(self, TimeLength):
		if len(self.l) < self.delta:
			self.E.append(self.E[self.k] + 0 - (TimeLength * self.Rs)) #TimeLength in h o s?
		else:
			self.E.append(self.E[self.k] + self.l[-self.delta] - (TimeLength * self.Rs)) #TimeLength in h o s?

	def computeDsBig(self, TimeLength):
		app = 0
		if len(self.l) < self.delta:
			app = (0 + self.E[self.k]) / TimeLength
			
		else:
			print("self.l[self.k]" + str(self.l[self.k]))
			print("self.E[self.k]" + str(self.E[self.k]))
			app = (self.l[-self.delta] + self.E[self.k]) / TimeLength #TimeLength in h o s?
			
		print("app: "+str(app))
		if(app > self.r_s_max):
			self.d_s_big = self.r_s_max
		else:
			self.d_s_big = app

	def computeL(self, TimeLength):
		self.l.append(self.l[self.k] + TimeLength * (self.Ss-self.Rs))

	def updateK(self, kappa):
		self.k=kappa
		

