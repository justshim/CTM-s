import cell as c
import supervisor as s
import factory as f
import xlrd 
import xlwt
import matplotlib.pyplot as plt
import numpy as np

#######################
# File read section:  #
#######################

## read file CTM_data from xls file

#loc = ("C:/A_Tesi/Python/CTM-s/CTM_data.xls")
loc = ("C:/A_Tesi/CTMs-identification/CTM_param_out.xls")
#loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")

wb = xlrd.open_workbook(loc)
sh = wb.sheet_by_index(0)
sh.cell_value(0,0)

wbt = xlwt.Workbook()
ws = wbt.add_sheet('Sheet1')

## read phi first cell from xls file

#loc_phi = ("C:/A_Tesi/Python/CTM-s/phi_1.xls")
loc_phi = ("C:/A_Tesi/CTMs-identification/CTM_param_out.xls")
#loc_phi = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/phi_1.xls")

wb_phi = xlrd.open_workbook(loc_phi)

#sh_phi = wb_phi.sheet_by_index(0) 		#sheet 0 is a "realistic" input (24h)
#sh_phi = wb_phi.sheet_by_index(1)		#sheet 1 is a "synthetic" input with 2 equal peaks (24h)
#sh_phi = wb_phi.sheet_by_index(2)		#sheet 2 is a "synthetic" input with 1 peak (24h)
#sh_phi = wb_phi.sheet_by_index(3)		#sheet 3 is a "synthetic" input with 1 peak (3h)
#sh_phi = wb_phi.sheet_by_index(4)		#sheet 4 is a flat input (24h)
sh_phi = wb_phi.sheet_by_index(2)

sh_phi.cell_value(0,0)

phi_zero=[]
for i in range(0, sh_phi.nrows):
	phi_zero.append(sh_phi.cell_value(i,0))

cell_in = 0
row = 0
duration = 8640 # k=24h=8640 , k=1h=360, k=3h=1080

####################################
# First iteration without stations #
####################################
fac = f.Factory()

## create the stretch via the factory
	# timeLength [h],   lastPhi,  phi_zero
fac.createStretch(10/3600, 5000, phi_zero) 

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

ttt = fac.stretches[0].computeTTT()

d0 = []

##################################
# Exectution of the simulation:  #
##################################

k=0
while k<duration: 
	
	fac.stretches[0].update(k)
	k = k + 1
	
	d0.append(fac.stretches[0].delta_big[k-1])

delta_max0 = max(d0)
print("max(delta_0): " + str(delta_max0))
print()

ws.write(0, 3, delta_max0)



for cell_in in range(1, sh.nrows-3, 1):
	cell_out = cell_in + 2
	for delta in range(30, 361, 30): #360k = 3600s
		
		###################################################
		# Initialization of all components of the model:  #
		###################################################

		## create a factory instance that manages the creation of objects
		fac = f.Factory()

		## create the stretch via the factory
			# timeLength [h],   lastPhi,  phi_zero
		fac.createStretch(10/3600, 5000, phi_zero) 

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

		#ttt = fac.stretches[0].computeTTT()


		## support variables to save various parameters during execution, and possibly plot them
		l0 = []
		e0 = []

		r0 = []
		r1 = []
		r2 = []
		r3 = []
		r4 = []
		r5 = []
		r6 = []
		r7 = []
		r8 = []
		r9 = []
		r10 = []
		r11 = []

		cong0 = []
		cong1 = []
		cong2 = []
		cong3 = []
		cong4 = []
		cong5 = []
		cong6 = []
		cong7 = []
		cong8 = []

		d = []

		##################################
		# Exectution of the simulation:  #
		##################################

		k=0
		while k<duration: 	
			#print("\nTime instant: " + str(k))
			
			fac.stretches[0].update(k)
			k = k + 1
			
			# save the various parameters in the previously created variables
			#l0.append(fac.stretches[0].stations[0].l[k])
			#e0.append(fac.stretches[0].stations[0].e[k])

			r0.append(fac.stretches[0].cells[0].rho[k])
			r1.append(fac.stretches[0].cells[1].rho[k])
			r2.append(fac.stretches[0].cells[2].rho[k])
			r3.append(fac.stretches[0].cells[3].rho[k])
			r4.append(fac.stretches[0].cells[4].rho[k])
			r5.append(fac.stretches[0].cells[5].rho[k])
			r6.append(fac.stretches[0].cells[6].rho[k])
			r7.append(fac.stretches[0].cells[7].rho[k])
			r8.append(fac.stretches[0].cells[8].rho[k])

			cong0.append(fac.stretches[0].cells[0].congestion_state)
			cong1.append(fac.stretches[0].cells[1].congestion_state)
			cong2.append(fac.stretches[0].cells[2].congestion_state)
			cong3.append(fac.stretches[0].cells[3].congestion_state)
			cong4.append(fac.stretches[0].cells[4].congestion_state)
			cong5.append(fac.stretches[0].cells[5].congestion_state)
			cong6.append(fac.stretches[0].cells[6].congestion_state)
			cong7.append(fac.stretches[0].cells[7].congestion_state)
			cong8.append(fac.stretches[0].cells[8].congestion_state)

			d.append(fac.stretches[0].delta_big[k-1])

		integ = 0
		for kkk in d:
			integ = integ + kkk
		
		print("i: " + str(cell_in) + "; j: " + str(cell_out))
		print("delta stazione :" + str(delta))
		print("Integ :" + str(integ))
		delta_max = max(d)
		pi = (delta_max0-delta_max)/delta_max0
		print("Pi :" + str(pi))
		print()

		row = row + 1
		ws.write(row, 0, cell_in)
		ws.write(row, 1, delta)
		ws.write(row, 2, integ)
		ws.write(row, 3, delta_max)
		ws.write(row, 4, pi)
		
		plt.figure(row)
		plt.grid(True)
		plt.plot(cong8)
		
	plt.show()	

#wbt.save('C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/opti_data.xls')	
wbt.save('C:/A_Tesi/Python/CTM-s/opti_data.xls')	
#print("Len rho: " + str(len(fac.stretches[0].cells[0].rho))) 

#############################
# Plot management section:  #
#############################

# plt.figure(0)
# plt.grid(True)
# plt.xlabel('k')
# plt.ylabel('r6')
# plt.plot(d)


