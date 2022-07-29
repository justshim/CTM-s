#!/usr/bin/env python
import cell as c
import supervisor as s
import factory as f
import xlrd 
import xlwt
import matplotlib.pyplot as plt
import numpy as np
import itertools
import multiprocessing

def mia_funz(params):
	
	global duration
	global delta_max0
	global timeLength
	global last_phi
	global phi_zero
	global length 
	global v_free 
	global w
	global q_max
	global rho_max
	###################################################
	# Initialization of all components of the model:  #
	###################################################

	## create a factory instance that manages the creation of objects
	fac = f.Factory()

	## create the stretch via the factory
		# timeLength [h],   lastPhi,  phi_zero
	fac.createStretch(timeLength, last_phi, phi_zero) 

	## create the cells via the factory
	for i in range(0, len(length)):
		           #ID stretch, length, v_free, w, q_max, rho_max, p_ms
		fac.addCellToStretch(0, length[i], v_free[i], w[i], q_max[i], rho_max[i], 0.95)	
				
	## create the stations via the factory
				#ID stretch, r_s_max, i, j, delta, beta_s, p
	fac.addStationToStretch(0, 500, params[0], params[1], params[2], 0.07, 0.05) 

	## create the on-ramps via the factory
				#ID stretch, d_r, r_r_max, j, p_r
	#fac.addOnRampToStretch(0, 100, 500, 2, 0.05)

	## create the off-ramps via the factory
				#ID_stretch, i, beta_r
	#fac.addOffRampToStretch(0, 7, 0.05)

	d=[]
	k=0
	while k<duration: #execution of the simulation with the current station settings
		fac.stretches[0].update(k)
		k = k + 1
		d.append(fac.stretches[0].delta_big[k-1])

	integ = 0
	for kkk in d:
		integ = integ + kkk
	
	delta_max = max(d)
	pi = (delta_max0-delta_max)/delta_max0
	tupla=[params[0], params[1], params[2], integ, delta_max, pi]
	return tupla

def main():
	a = range(1, sh.nrows-1, 1)
	b = range(1, sh.nrows, 1)
	c = range(60,721, 60)
	d = range(10)
	paramlist = list(itertools.product(a,b,c)) #all possible combinations of tuples
	
	paramlist2 = []
	for i in range(len(paramlist)):
		if(paramlist[i][1] > paramlist[i][0]): #only if exit station > enter station
			paramlist2.append(paramlist[i])

	pool = multiprocessing.Pool() #generate processes equal to the number of cores
	res = pool.map(mia_funz,paramlist2) #distribute the parameter sets evenly across the cores
	return res 
	
#######################
# File read section:  #
#######################
#path_file_output = 'C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/opti_data.xls'
#path_file_output = 'C:/A_Tesi/Python/CTM-s/opti_data.xls'
path_file_output = ("H:/Il mio Drive/Tesi magistrale/Python/CTM-s/opti_data.xls")
wbt = xlwt.Workbook()
out_file =  wbt.add_sheet('Sheet1')
## read file CTM_data from xls file

#loc = ("C:/A_Tesi/Python/CTM-s/CTM_data.xls")
#loc = ("C:/A_Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
loc = ("H:/Il mio Drive/Tesi magistrale/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
#loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")

input_file = xlrd.open_workbook(loc)
sh = input_file.sheet_by_name("Cells parameters")
sh.cell_value(0,0)

sh_phi = input_file.sheet_by_name("First Demand Smooth")
sh_phi.cell_value(0,0)
phi_zero=[]
for i in range(0, sh_phi.nrows):
	phi_zero.append(sh_phi.cell_value(i,0))

sh_last_phi = input_file.sheet_by_name("Last Demand Smooth")
sh_last_phi.cell_value(0,0)
last_phi=[]
for i in range(0, sh_last_phi.nrows):
	last_phi.append(sh_last_phi.cell_value(i,0))

duration = 8640 # k=24h=8640 , k=1h=360, k=3h=1080
length=[]
v_free=[]
w=[]
q_max=[]
rho_max=[]
timeLength = sh.cell_value(2,6)
for i in range(1, sh.nrows):
	length.append(sh.cell_value(i,1))
	v_free.append(sh.cell_value(i,2))
	w.append(sh.cell_value(i,3))
	q_max.append(sh.cell_value(i,4))
	rho_max.append(sh.cell_value(i,5))

####################################
# First iteration without stations #
####################################
fac = f.Factory()

## create the stretch via the factory
	# timeLength [h],   lastPhi,  phi_zero
fac.createStretch(timeLength, last_phi, phi_zero) 

## create the cells via the factory
for i in range(0, len(length)):
		           #ID stretch, length, v_free, w, q_max, rho_max, p_ms
		fac.addCellToStretch(0, length[i], v_free[i], w[i], q_max[i], rho_max[i], 0.95)	
## create the on-ramps via the factory
			#ID stretch, d_r, r_r_max, j, p_r
#fac.addOnRampToStretch(0, 100, 500, 2, 0.05)

## create the off-ramps via the factory
			#ID_stretch, i, beta_r
#fac.addOffRampToStretch(0, 7, 0.05)

ttt = fac.stretches[0].computeTTT() #compute TTT
d0=[] 
k=0
while k<duration: #execution of the first simulation without any station
	fac.stretches[0].update(k)
	k = k + 1
	d0.append(fac.stretches[0].delta_big[k-1])

delta_max0 = max(d0)
if __name__ == "__main__":    
	res=main()
	
	#Writing output file
	out_file.write(0, 0, "i")
	out_file.write(0, 1, "j")
	out_file.write(0, 2, "delta")
	#out_file.write(0, 3, "beta")
	out_file.write(0, 3, "integral")
	out_file.write(0, 4, "max_delta")
	out_file.write(0, 5, "pi")

	for i in range(1, len(res)):
		out_file.write(i, 0, res[i][0])
		out_file.write(i, 1, res[i][1])
		out_file.write(i, 2, res[i][2])
		out_file.write(i, 3, res[i][3])
		out_file.write(i, 4, res[i][4])
		out_file.write(i, 5, res[i][5])

	wbt.save(path_file_output)