from typing import List
from model.supervisor import Stretch
from model.parameters import CTMsParameters
import numpy as np


class Factory:
	"""
	Class for the creation of instances of cells and stations
	"""

	stretches: List[Stretch]
	n_stretches: int

	def __init__(self, parameters: CTMsParameters, data: np.ndarray):
		self.stretches = []
		self.n_stretches = 0
		self.build(parameters, data)

	def create_stretch(self, dt, first_phi, last_phi):
		"""
		Method for the creation of instances of the object Stretch
		"""

		stretch = Stretch(dt, first_phi, last_phi)
		self.stretches.append(stretch)
		self.n_stretches += 1

	def add_cell_to_stretch(self, id_stretch, length, v_free, w, q_max, rho_max, p):
		"""
		Call to supervisors' method for the creation of instances of the object Cell
		"""

		self.stretches[id_stretch].createCell(length, v_free, w, q_max, rho_max, p)

	def add_station_to_stretch(self, id_stretch, r_s_max, i, j, delta, beta_s, p):
		"""
		Call to supervisors' method for the creation of instances of the object Station
		"""

		self.stretches[id_stretch].createStation(r_s_max, i, j, delta, beta_s, p)

	def add_on_ramp_to_stretch(self, id_stretch, d_r, r_r_max, j, p_r):
		"""
		Call to supervisors' method for the creation of instances of the object OnRamp
		"""

		self.stretches[id_stretch].createOnRamp(d_r, r_r_max, j, p_r)

	def add_off_ramp_to_stretch(self, id_stretch, i, beta_r):
		"""
		Call to supervisors' method for the creation of instances of the object OffRamp
		"""

		self.stretches[id_stretch].createOffRamp(i, beta_r)

	def build(self, parameters: CTMsParameters, phi_data: np.ndarray):
		"""
		Initialization helper function to build factory according to parameters
		"""
		# Create highway stretch
		self.create_stretch(
			dt=parameters.highway.dt[0],
			first_phi=phi_data,
			last_phi=phi_data)

		# Add cells to highway stretch
		for i in range(len(parameters.highway)):
			self.add_cell_to_stretch(
				id_stretch=0,
				length=parameters.highway.l[i],
				v_free=parameters.highway.v[i],
				w=parameters.highway.w[i],
				q_max=parameters.highway.q_max[i],
				rho_max=parameters.highway.rho_max[i],
				p=parameters.highway.p_ms[i])

		# Add on-ramps to the stretch
		for i in range(len(parameters.onramps)):
			self.add_on_ramp_to_stretch(
				id_stretch=0,
				d_r=parameters.onramps.d_r,
				r_r_max=parameters.onramps.r_r_max[i],
				j=parameters.onramps.j[i],
				p_r=parameters.onramps.p_r[i])

		# Add off-ramps to the stretch
		for i in range(len(parameters.offramps)):
			self.add_off_ramp_to_stretch(
				id_stretch=0,
				i=parameters.offramps.i[i],
				beta_r=parameters.offramps.beta_r[i])

		# Add service station to the stretch
		for i in range(len(parameters.stations)):
			self.add_station_to_stretch(
				id_stretch=0,
				r_s_max=parameters.stations.r_s_max[i],
				i=parameters.stations.i[i],
				j=parameters.stations.j[i],
				delta=parameters.stations.delta[i],
				beta_s=parameters.stations.beta_s[i],
				p=parameters.stations.p[i])
