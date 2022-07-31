#!/usr/bin/env python

from model import cell as c
from model import supervisor as s
from model import factory as f
import xlrd 
import xlwt
import matplotlib.pyplot as plt
import numpy as np
import itertools
import multiprocessing
import os
import csv
from tqdm.auto import tqdm
from p_tqdm import p_map, p_umap, p_imap, p_uimap

class data_creation:
	
	def __init__(self, loc, path_file_output, duration):		
		self.path_file_output = path_file_output
		self.loc = loc
		self.phi_zero=[]
		self.last_phi=[]
		self.duration = duration 
		self.length=[]
		self.v_free=[]
		self.w=[]
		self.q_max=[]
		self.rho_max=[]
		self.timeLength = 0
		self.delta_max0 = 0
		self.nrow = 0

	def initVars(self):
		input_file = xlrd.open_workbook(self.loc)

		sh = input_file.sheet_by_name("Cells parameters")
		sh.cell_value(0,0)
		self.nrow = sh.nrows

		sh_phi = input_file.sheet_by_name("First Demand Smooth")
		sh_phi.cell_value(0,0)
		
		sh_last_phi = input_file.sheet_by_name("Last Demand Smooth")
		sh_last_phi.cell_value(0,0)
		
		for i in range(0, sh_phi.nrows):
			self.phi_zero.append(sh_phi.cell_value(i,0))

		for i in range(0, sh_last_phi.nrows):
			self.last_phi.append(sh_last_phi.cell_value(i,0))
		
		self.timeLength = sh.cell_value(2,6)
		for i in range(1, sh.nrows):
			self.length.append(sh.cell_value(i,1))
			self.v_free.append(sh.cell_value(i,2))
			self.w.append(sh.cell_value(i,3))
			self.q_max.append(sh.cell_value(i,4))
			self.rho_max.append(sh.cell_value(i,5))

		pass

	def firstIteration(self):
		fac = f.Factory()

		## create the stretch via the factory
			# timeLength [h],   lastPhi,  phi_zero
		fac.createStretch(self.timeLength, self.last_phi, self.phi_zero) 

		## create the cells via the factory
		for i in range(0, len(self.length)):
				           #ID stretch, length, v_free, w, q_max, rho_max, p_ms
				fac.addCellToStretch(0, self.length[i], self.v_free[i], self.w[i], self.q_max[i], self.rho_max[i], 1)	
		## create the on-ramps via the factory
					#ID stretch, d_r, r_r_max, j, p_r
		#fac.addOnRampToStretch(0, 100, 500, 2, 0.05)

		## create the off-ramps via the factory
					#ID_stretch, i, beta_r
		#fac.addOffRampToStretch(0, 7, 0.05)

		ttt = fac.stretches[0].computeTTT() #compute TTT
		d0=[] 
		k=0
		while k<self.duration: #execution of the first simulation without any station
			fac.stretches[0].update(k)
			k = k + 1
			d0.append(fac.stretches[0].delta_big[k-1])

		self.delta_max0 = max(d0)
		pass
	
	def initSimulation(self, initDelta, delta, stepDelta, initBeta, beta, stepBeta, initPriority, priority, stepPriority):
		a = range(1, self.nrow-1, 1)
		b = range(1, self.nrow, 1)
		c = range(initDelta,delta, stepDelta)
		d = range(initBeta, beta, stepBeta)
		e = range(initPriority,priority,stepPriority)

		paramlist = list(itertools.product(a,b,c,d,e)) #all possible combinations of tuples
		
		paramlist2 = []
		for i in range(len(paramlist)):
			if(paramlist[i][1] > paramlist[i][0]): #only if exit station > enter station
				paramlist2.append(paramlist[i])
		print("Total cases to be evaluated: " + str(len(paramlist2)))
		
		with multiprocessing.Pool(processes=os.cpu_count()) as pool:
		#pool = multiprocessing.Pool() #generate processes equal to the number of cores
			chuncksize=int(len(paramlist2)/os.cpu_count())
			res = p_map(self.simulation,paramlist2)

		
		#res = pool.map(self.simulation,paramlist2) #distribute the parameter sets evenly across the cores
		#pool.close() # No more work
		#pool.join() # Wait for completion
		self.writeResults(res)
		pass 

	def simulation(self, params):
		beta=params[3]/100
		p_ms=params[4]/100
		p_station=(100-params[4])/100
		
		## create a factory instance that manages the creation of objects
		fac = f.Factory()

		## create the stretch via the factory
							#timeLength [h],   lastPhi,  phi_zero
		fac.createStretch(self.timeLength, self.last_phi, self.phi_zero) 

		## create the cells via the factory
		for i in range(0, len(self.length)):
			           #ID stretch, length, 	v_free, 			w, 			q_max, 			rho_max, 	p_ms
			fac.addCellToStretch(0, self.length[i], self.v_free[i], self.w[i], self.q_max[i], self.rho_max[i], p_ms)	
					
		## create the stations via the factory
					#ID stretch, r_s_max, i, j, delta, beta_s, p
		fac.addStationToStretch(0, 500, params[0], params[1], params[2], beta, p_station) 

		## create the on-ramps via the factory
					#ID stretch, d_r, r_r_max, j, p_r
		#fac.addOnRampToStretch(0, 100, 500, 2, 0.05)

		## create the off-ramps via the factory
					#ID_stretch, i, beta_r
		#fac.addOffRampToStretch(0, 7, 0.05)

		d=[]
		k=0
		while k<self.duration: #execution of the simulation with the current station settings
			fac.stretches[0].update(k)
			k = k + 1
			d.append(fac.stretches[0].delta_big[k-1])

		integ = 0
		for kkk in d:
			integ = integ + kkk
		
		delta_max = max(d)
		pi = (self.delta_max0-delta_max)/self.delta_max0
		tupla=[params[0], params[1], params[2], beta, p_station, integ, delta_max, pi]
		return tupla
	
	def writeResults(self, res):
		print("Writing output file... ")
		
		with open(self.path_file_output, 'w', encoding='UTF8', newline='') as f:
			writer = csv.writer(f)
			header =["i","j","delta","beta","priority","integral","max_delta","pi"]
			writer.writerow(header)
			for i in range(0, len(res)):
				writer.writerow(res[i])
		print("Done!")
		pass