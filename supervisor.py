import cell as c
import station as st

class Stretch:
	"""Controller class for the model, represents the system at large"""

	def __init__(self):
		self.cells = []
		self.stations = []
		self.n_cells = 0
		self.n_stations = 0

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

	def createStation(self, r_s_max, i, j):
		station = st.Station(self.n_stations, r_s_max, i, j) 
		self.stations.append(station)
		self.n_stations = self.n_stations + 1