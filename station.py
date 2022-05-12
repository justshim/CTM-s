import service as ser
import supervisor as su

class Station:
	"""Class modeling the service stations on an highway stretch in the CTM-s model"""
	def __init__(self, ID, r_s_max, i, j):
		self.ID = ID
		self.r_s_max = r_s_max
		self.i = i
		self.j = j
		self.Ss = 0
		self.Rs = 0
		self.E = 0
		self.d_s_big = 0
		self.services = []
		self.n_services = 0

	def toString(self):
		print("Station ID: "+str(self.ID))
		print("From cell "+str(self.i))
		print("To cell "+str(self.j))

	def computeSs(self, phiMeno):
		total_beta=0
		for serv in self.services:
			total_beta=total_beta+serv.beta_s
		return total_beta*phiMeno

	def computeRs(self):
		pass

	def computeE(self, TimeLength):
		total_l=0
		for serv in self.services:
			total_l=total_l+serv.l
		self.E=self.E+total_l-TimeLength*self.Rs

	def computeDsBig(self):
		#DA VERIFICARE SE RS è CORRETTO
		if(self.Rs > self.r_s_max):
			self.d_s_big=self.r_s_max
		else:
			self.d_s_big=self.Rs

	def createService(self, delta, beta_s):
		service = ser.Service(self.n_services, delta, beta_s) 
		self.services.append(service)
		self.n_services = self.n_services + 1

	def updateService(self, TimeLength):
		for serv in self.services:
			serv.computeL(TimeLength, self.Ss, self.Rs)

	def computeTotalBeta():
		beta=0
		for serv in self.services:
			beta=beta+serv.beta_s
		return beta
