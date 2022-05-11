import service as ser

class Station:
	"""Class modeling the service stations on an highway stretch in the CTM-s model"""
	def __init__(self, ID, r_s_max, i, j):
		self.ID = ID
		self.r_s_max = r_s_max
		self.i = i
		self.j = j

		self.services = []
		self.n_services = 0

	def toString(self):
		print("Station ID: "+str(self.ID))
		print("From cell "+str(self.i))
		print("To cell "+str(self.j))

	def computeSs(self):
		pass

	def computeRs(self):
		pass

	def computeE(self):
		pass

	def computeDsBig(self):
		pass

	def createService(self, delta, beta_s):
		service = ser.Service(self.n_services, delta, beta_s) 
		self.services.append(service)
		self.n_services = self.n_services + 1