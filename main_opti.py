import cell as c
import supervisor as s
import factory as f
import xlrd 
import xlwt
import matplotlib.pyplot as plt
import numpy as np



def mia_funz(sh, duration, last_phi, phi_zero, delta_max0):
	row=0
	bbb=[]
	#setup any possible station with any possible settings
	for cell_in in range(1, sh.nrows-1, 1):

		for cell_out in range(cell_in+1, sh.nrows, 1):

			for delta in range(60,61, 60): #360k = 3600s = 60 min. 720k = 7200s = 120 min
				
				#for beta in range(1, 21, 1):
					#beta=beta/100
					###################################################
					# Initialization of all components of the model:  #
					###################################################

					## create a factory instance that manages the creation of objects
					fac = f.Factory()

					## create the stretch via the factory
						# timeLength [h],   lastPhi,  phi_zero
					fac.createStretch(sh.cell_value(2,6), last_phi, phi_zero) 

					## create the cells via the factory
					for i in range(1, sh.nrows):
						           #ID stretch   length,             v_free,                   w,              ,   q_max,           rho_max,         p_ms
						fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), sh.cell_value(i,4), sh.cell_value(i,5), 0.95)	
								
					## create the stations via the factory
								#ID stretch, r_s_max, i, j, delta, beta_s, p
					fac.addStationToStretch(0, 500, cell_in, cell_out, delta, 0.07, 0.05) 

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

					row = row + 1
					print(str(row)+"/936")
					aaa=[row, cell_in, cell_out, delta, integ, delta_max, pi]
					bbb.append(aaa)
	return bbb



def main():
	    
	#######################
	# File read section:  #
	#######################
	#path_file_output = 'C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/opti_data.xls'
	#path_file_output = 'C:/A_Tesi/Python/CTM-s/opti_data.xls'
	path_file_output = ("H:/Il mio Drive/Tesi magistrale/Python/CTM-s/opti_data.xls")

	## read file CTM_data from xls file

	#loc = ("C:/A_Tesi/Python/CTM-s/CTM_data.xls")
	#loc = ("C:/A_Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
	loc = ("H:/Il mio Drive/Tesi magistrale/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
	#loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")

	wb = xlrd.open_workbook(loc)
	sh = wb.sheet_by_name("Cells parameters")
	sh.cell_value(0,0)

	wbt = xlwt.Workbook()
	ws = wbt.add_sheet('Sheet1')

	## read phi first cell from xls file

	#loc_phi = ("C:/A_Tesi/Python/CTM-s/phi_1.xls")
	#loc_phi = ("C:/A_Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
	#loc_phi = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/phi_1.xls")
	loc_phi = loc
	wb_phi = xlrd.open_workbook(loc_phi)

	#sh_phi = wb_phi.sheet_by_index(0) 		#sheet 0 is a "realistic" input (24h)
	#sh_phi = wb_phi.sheet_by_index(1)		#sheet 1 is a "synthetic" input with 2 equal peaks (24h)
	#sh_phi = wb_phi.sheet_by_index(2)		#sheet 2 is a "synthetic" input with 1 peak (24h)
	#sh_phi = wb_phi.sheet_by_index(3)		#sheet 3 is a "synthetic" input with 1 peak (3h)
	#sh_phi = wb_phi.sheet_by_index(4)		#sheet 4 is a flat input (24h)

	sh_phi = wb_phi.sheet_by_name("First Demand Smooth")
	sh_phi.cell_value(0,0)
	phi_zero=[]
	for i in range(0, sh_phi.nrows):
		phi_zero.append(sh_phi.cell_value(i,0))

	sh_last_phi = wb.sheet_by_name("Last Demand Smooth")
	sh_last_phi.cell_value(0,0)
	last_phi=[]
	for i in range(0, sh_last_phi.nrows):
		last_phi.append(sh_last_phi.cell_value(i,0))


	duration = 8640 # k=24h=8640 , k=1h=360, k=3h=1080

	####################################
	# First iteration without stations #
	####################################
	fac = f.Factory()

	## create the stretch via the factory
		# timeLength [h],   lastPhi,  phi_zero
	fac.createStretch(sh.cell_value(2,6), last_phi, phi_zero) 

	## create the cells via the factory
	for i in range(1, sh.nrows):
		           #ID stretch   length,             v_free,                   w,              ,   q_max,           rho_max,         p_ms
		fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), sh.cell_value(i,4), sh.cell_value(i,5), 0.95)	

	## create the on-ramps via the factory
				#ID stretch, d_r, r_r_max, j, p_r
	#fac.addOnRampToStretch(0, 100, 500, 2, 0.05)

	## create the off-ramps via the factory
				#ID_stretch, i, beta_r
	#fac.addOffRampToStretch(0, 7, 0.05)

	ttt = fac.stretches[0].computeTTT() #compute TTT before any simulation for reference
	d0=[] 
	k=0
	while k<duration: #execution of the first simulation without any station
		fac.stretches[0].update(k)
		k = k + 1
		d0.append(fac.stretches[0].delta_big[k-1])

	delta_max0 = max(d0)
	#print("max(delta_0): " + str(delta_max0))
	#print()


	bb=mia_funz(sh, duration, last_phi, phi_zero, delta_max0)

	#Writing output file
	ws.write(0, 0, "i")
	ws.write(0, 1, "j")
	ws.write(0, 2, "delta")
	#ws.write(0, 3, "beta")
	ws.write(0, 3, "integral")
	ws.write(0, 4, "max_delta")
	ws.write(0, 5, "pi")

	for i in range(0, len(bb)):
		ws.write(bb[i][0], 0, bb[i][1])
		ws.write(bb[i][0], 1, bb[i][2])
		ws.write(bb[i][0], 2, bb[i][3])
		ws.write(bb[i][0], 3, bb[i][4])
		ws.write(bb[i][0], 4, bb[i][5])
		ws.write(bb[i][0], 5, bb[i][6])

	wbt.save(path_file_output)	


if __name__ == "__main__":
    main()