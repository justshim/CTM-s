class Service:
	"""Class modeling the services in a service station on an highway stretch in the CTM-s model"""
	def __init__(self, ID, delta, beta_s):
		self.ID = ID
		self.delta = delta
		self.beta_s = beta_s


	def toString(self):
		print("Service ID: "+str(self.ID))
		print("Time delay: "+str(self.delta))
		print("Split ratio: "+str(self.beta_s))

	def computeL():
		pass
		
		