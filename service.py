class Service:
	"""Class modeling the services in a service station on an highway stretch in the CTM-s model"""
	def __init__(self, ID, delta, beta_s):
		self.ID = ID
		self.delta = delta
		self.beta_s = beta_s
		self.l = 0

	def toString(self):
		print("Service ID: "+str(self.ID))
		print("Time delay: "+str(self.delta))
		print("Split ratio: "+str(self.beta_s))

	def computeL(self, TimeLength, Ss, Rs):
		# ATTENZIONE DA RIVEDERE: Ss e Rs sono riferiti alla stazione e non al servizio
		self.l=self.l+ TimeLength*(Ss*self.beta_s-Rs*self.beta_s)
		
		