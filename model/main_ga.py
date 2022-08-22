from model import cell as c
from model import supervisor as s
from model import factory as f
import xlrd 
import matplotlib.pyplot as plt
import numpy as np


#############################################
# input = path di ctm-param-out-niceee
# beta, priority, delta, i e j come variabili che variano in automatico dal G.A.
# output = integral delta e pi
#############################################


#######################
# File read section:  #
#######################

## read file CTM_data from xls file

loc = ("H:/Il mio Drive/Tesi magistrale/CTMs-identification/fnc/extracted_data/CTM_param_out_nice.xls")
#loc = ("C:/A_Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
#loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
#loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")

wb = xlrd.open_workbook(loc)
sh = wb.sheet_by_name("Cells parameters")
sh.cell_value(0,0)

## read phi first cell from xls file

#loc_phi = ("C:/A_Tesi/Python/CTM-s/phi_1.xls")
#loc_phi = ("C:/A_Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
#loc_phi = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
#loc_phi = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/phi_1.xls")
loc_phi=loc
wb_phi = xlrd.open_workbook(loc_phi)

#sh_phi = wb_phi.sheet_by_index(0) 		#sheet 0 is a "realistic" input (24h)
#sh_phi = wb_phi.sheet_by_index(1)		#sheet 1 is a "synthetic" input with 2 equal peaks (24h)
#sh_phi = wb_phi.sheet_by_index(2)		#sheet 2 is a "synthetic" input with 1 peak (24h)
#sh_phi = wb_phi.sheet_by_index(3)		#sheet 3 is a "synthetic" input with 1 peak (3h)
#sh_phi = wb_phi.sheet_by_index(4)		#sheet 4 is a flat input (24h)
#sh_phi = wb_phi.sheet_by_index(5)		#sheet 5 is a real input with peak at 1800 veh/h (24h)
#sh_phi = wb_phi.sheet_by_index(6)		#sheet 6 is a real input with peak at 2500 veh/h (24h)

sh_phi = wb_phi.sheet_by_name("First Demand Real")
sh_phi.cell_value(0,0)
phi_zero=[]
for i in range(0, sh_phi.nrows):
	phi_zero.append(sh_phi.cell_value(i,0))

sh_last_phi = wb.sheet_by_name("Last Demand Real")
sh_last_phi.cell_value(0,0)
last_phi=[]
for i in range(0, sh_last_phi.nrows):
	last_phi.append(sh_last_phi.cell_value(i,0))

###################################################
# Initialization of all components of the model:  #
###################################################

## create a factory instance that manages the creation of objects
fac = f.Factory()

## create the stretch via the factory
		 		#timeLength [h],   lastPhi, 	phi_zero
fac.createStretch(sh.cell_value(2,6), last_phi, phi_zero) 

## create the cells via the factory
for i in range(1, sh.nrows):
	           #ID stretch   length,             v_free,                   w,              ,   q_max,           rho_max,         p_ms
	fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), sh.cell_value(i,4), sh.cell_value(i,5), 1)	
		
## create the stations via the factory
			#ID stretch, r_s_max, i, j, delta, beta_s, p
#fac.addStationToStretch(0, 500, 3, 6, 60, 0.07, 0.05) #Note: r_s_max was statically assigned to the Qmax(4)/10 (the cell where the station merges back)

## create the on-ramps via the factory
			#ID stretch, d_r, r_r_max, j, p_r
#fac.addOnRampToStretch(0, 100, 500, 2, 0.05)

## create the off-ramps via the factory
			#ID_stretch, i, beta_r
#fac.addOffRampToStretch(0, 7, 0.05)

fac.stretches[0].computeTTT()

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
r12 = []

s0 = []
s1 = []
s2 = []
s3 = []
s4 = []
s5 = []
s6 = []
s7 = []
s8 = []
s9 = []
s10 = []
s11 = []
s12 = []

v10 = []

cong0 = []
cong1 = []
cong2 = []
cong3 = []
cong4 = []
cong5 = []
cong6 = []
cong7 = []
cong8 = []
cong9 = []
cong10 = []
cong11 = []

d = []
demand = []
supply = []
demand_next = []
supply_next = []
fipiu = []
fimeno = []
##################################
# Exectution of the simulation:  #
##################################

k=0
while k<8640: 	# k=24h=8640 , k=1h=360, k=3h=1080
	#print("\nTime instant: " + str(k))
	
	fac.stretches[0].update(k)
	k = k + 1

	d.append(fac.stretches[0].delta_big[k-1])

print("\nEnd")
	
#print("Len rho: " + str(len(fac.stretches[0].cells[0].rho))) 

#############################
# Plot management section:  #
#############################

plt.figure(0)
plt.grid(True)
plt.xlabel('k')
plt.ylabel('r0')
plt.plot(r0)

plt.figure(1)
plt.grid(True)
plt.xlabel('k')
plt.ylabel('rho')
plt.plot(r9)
plt.plot(r10)

plt.figure(2)
plt.grid(True)
plt.xlabel('k')
plt.ylabel('v10')
plt.plot(v10)

plt.figure(3)
plt.grid(True)
plt.xlabel('k')
plt.ylabel('r12')
plt.plot(r12)

# plt.figure(4)
# plt.grid(True)
# plt.xlabel('k')
# plt.ylabel('cong0')
# plt.plot(cong0)

# #plt.figure(5)
# plt.grid(True)
# plt.xlabel('k')
# plt.ylabel('cong1')
# plt.plot(cong1)

plt.figure(6)
plt.grid(True)
plt.xlabel('k')
plt.ylabel('congestion')
plt.plot(cong9)
plt.plot(cong10)

plt.figure(7)
plt.grid(True)
plt.xlabel('k')
plt.plot(demand)
plt.plot(demand_next)
plt.legend(['demand','demand_next'], title = "Legend")

plt.figure(77)
plt.grid(True)
plt.xlabel('k')
plt.plot(s0)
plt.plot(s1)
plt.plot(s2)
plt.plot(s3)
plt.plot(s4)
plt.plot(s5)
plt.plot(s6)
plt.plot(s7)
plt.plot(s8)
plt.plot(s9)
plt.plot(s10)
plt.plot(s11)
plt.plot(s12)
plt.legend(['s0','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s11','s12'], title = "Legend")

plt.figure(8)
plt.grid(True)
plt.xlabel('k')
plt.ylabel('fipiu')
plt.plot(fipiu)


plt.grid(True)
plt.xlabel('k')
plt.ylabel('fimeno')
plt.plot(fimeno)


# plt.figure(7)
# plt.grid(True)
# plt.xlabel('k')
# plt.ylabel('cong10')
# plt.plot(cong10)

# plt.figure(99)
# plt.grid(True)
# plt.xlabel('k')
# plt.ylabel('cong11')
# plt.plot(cong11)

# plt.figure(99)
# plt.grid(True)
# plt.plot(phi_zero)
# plt.scatter(np.linspace(0, 8640, 8640), l0, 4, marker="x")

plt.show()