from model import supervisor as s


class Factory:
	# TODO: Provide type hints for clarity
	# TODO: Clean up method function names
	# TODO: Address difference in Python indexing and cell ID indices...

	"""
	Class for the creation of instances of cells and stations
	"""
	def __init__(self):
		self.stretches = []
		self.n_stretches = 0

	def createStretch(self, time_length, last_phi, first_d_big):
		"""
		Method for the creation of instances of the object Stretch
		"""

		stretch = s.Stretch(time_length, last_phi, first_d_big)
		self.stretches.append(stretch)
		self.n_stretches += 1

	def addCellToStretch(self, id_stretch, length, v_free, w, q_max, rho_max, p):
		"""
		Call to supervisors' method for the creation of instances of the object Cell
		"""

		self.stretches[id_stretch].createCell(length, v_free, w, q_max, rho_max, p)

	def addStationToStretch(self, id_stretch, r_s_max, i, j, delta, beta_s, p):
		"""
		Call to supervisors' method for the creation of instances of the object Station
		"""

		self.stretches[id_stretch].createStation(r_s_max, i, j, delta, beta_s, p)

	def addOnRampToStretch(self, id_stretch, d_r, r_r_max, j, p_r):
		"""
		Call to supervisors' method for the creation of instances of the object OnRamp
		"""

		self.stretches[id_stretch].createOnRamp(d_r, r_r_max, j, p_r)

	def addOffRampToStretch(self, id_stretch, i, beta_r):
		"""
		Call to supervisors' method for the creation of instances of the object OffRamp
		"""

		self.stretches[id_stretch].createOffRamp(i, beta_r)
		