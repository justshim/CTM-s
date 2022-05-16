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
		self.E = 0
		self.oldE = 0
		self.d_s_big = 0
		self.delta = delta
		self.beta_s = beta_s
		self.l = [0]
		self.p = p

	def toString(self):
		print("Station ID: "+str(self.ID_station))
		print("From cell "+str(self.i))
		print("To cell "+str(self.j))
		print("Time delay: "+str(self.delta))
		print("Split ratio: "+str(self.beta_s))
		print("Vehicles number: "+str(self.l))
		print()

	def computeSs(self, phiMeno):
		self.Ss = self.beta_s*phiMeno

	def computeRs(self):
			self.Rs=self.d_s_big
			print("d_s_big " + str(self.d_s_big))

	def computeE(self, TimeLength):
		self.oldE = self.E
		
		if len(self.l) < self.delta:
			self.E = self.oldE + 0 - (TimeLength * self.Rs)
		else:
			self.E = self.oldE + self.l[-self.delta] - (TimeLength * self.Rs)

	def computeDsBig(self, TimeLength):
		app = 0
		if len(self.l) < self.delta:
			app = (0 + self.oldE) / TimeLength
		else:
			app = (self.l[-self.delta] + self.oldE) / TimeLength

		if(app > self.r_s_max):
			self.d_s_big = self.r_s_max
		else:
			self.d_s_big = app

	def computeL(self, TimeLength):
		new_l = self.l[-1] + TimeLength * (self.Ss-self.Rs)
		self.l.append(new_l)

