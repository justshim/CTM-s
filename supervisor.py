import cell as c

class Stretch:
	"""Controller class for the model, represents the system at large"""
	
	def __init__(self):
		self.cells = []
		self.n_cells = 0

	def computeEcal(self):
		pass

	def computeDelta(self):
		pass

	def computePhi(self):
		pass
		
	def createCell(self, length, v, w, q, q_max, s, r, rho_max, beta, p):
		cell = c.Cell(self.n_cells, length, v, w, q, q_max, s, r, rho_max, beta, p) 
		self.cells.append(cell)
		self.n_cells = self.n_cells + 1