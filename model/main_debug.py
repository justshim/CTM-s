from model import cell as c
from model import supervisor as s
from model import factory as f
import xlrd 
import matplotlib.pyplot as plt
import numpy as np

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

# phi_1_24h_doublepeak
# phi_1_24h_real
# phi_1_24h_realsmooth
# phi_1_24h_realsmooth_incr
# phi_1_24h_singlepeak

path_phi = "C:\A_Tesi\Aimsun\CTM-s\phi_input\phi_1_24h_realsmooth_incr.txt"
phi_zero=[]
last_phi=[]
with open(path_phi, encoding='utf8') as fa:
	for line in fa:
		phi_zero.append(float(line.strip()))
		last_phi.append(float(line.strip()))


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
fac.addStationToStretch(0, 500, 8, 10, 240, 0.07, 0.05)

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
phi = []
##################################
# Exectution of the simulation:  #
##################################

k=0
while k<8640: 	# k=24h=8640 , k=1h=360, k=3h=1080
	#print("\nTime instant: " + str(k))
	
	fac.stretches[0].update(k)
	k = k + 1
	# save the various parameters in the previously created variables
	l0.append(fac.stretches[0].stations[0].l[k])
	e0.append(fac.stretches[0].stations[0].e[k])

	r0.append(fac.stretches[0].cells[0].rho[k])
	r1.append(fac.stretches[0].cells[1].rho[k])
	r2.append(fac.stretches[0].cells[2].rho[k])
	r3.append(fac.stretches[0].cells[3].rho[k])
	r4.append(fac.stretches[0].cells[4].rho[k])
	r5.append(fac.stretches[0].cells[5].rho[k])
	r6.append(fac.stretches[0].cells[6].rho[k])
	r7.append(fac.stretches[0].cells[7].rho[k])
	r8.append(fac.stretches[0].cells[8].rho[k])
	r9.append(fac.stretches[0].cells[9].rho[k])
	r10.append(fac.stretches[0].cells[10].rho[k])
	r11.append(fac.stretches[0].cells[11].rho[k])
	r12.append(fac.stretches[0].cells[12].rho[k])

	s0.append(fac.stretches[0].cells[0].s_big)
	s1.append(fac.stretches[0].cells[1].s_big)
	s2.append(fac.stretches[0].cells[2].s_big)
	s3.append(fac.stretches[0].cells[3].s_big)
	s4.append(fac.stretches[0].cells[4].s_big)
	s5.append(fac.stretches[0].cells[5].s_big)
	s6.append(fac.stretches[0].cells[6].s_big)
	s7.append(fac.stretches[0].cells[7].s_big)
	s8.append(fac.stretches[0].cells[8].s_big)
	s9.append(fac.stretches[0].cells[9].s_big)
	s10.append(fac.stretches[0].cells[10].s_big)
	s11.append(fac.stretches[0].cells[11].s_big)
	s12.append(fac.stretches[0].cells[12].s_big)

	phi.append(fac.stretches[0].cells[9].phi)

	demand.append(fac.stretches[0].cells[9].d_big)
	supply.append(fac.stretches[0].cells[9].s_big)
	demand_next.append(fac.stretches[0].cells[10].d_big)
	supply_next.append(fac.stretches[0].cells[10].s_big)

	fipiu.append(fac.stretches[0].cells[10].phi_plus)
	fimeno.append(fac.stretches[0].cells[10].phi_minus)

	v10.append(fac.stretches[0].cells[10].v[k-1])

	cong0.append(fac.stretches[0].cells[0].congestion_state)
	cong1.append(fac.stretches[0].cells[1].congestion_state)
	cong2.append(fac.stretches[0].cells[2].congestion_state)
	cong3.append(fac.stretches[0].cells[3].congestion_state)
	cong4.append(fac.stretches[0].cells[4].congestion_state)
	cong5.append(fac.stretches[0].cells[5].congestion_state)
	cong6.append(fac.stretches[0].cells[6].congestion_state)
	cong7.append(fac.stretches[0].cells[7].congestion_state)
	cong8.append(fac.stretches[0].cells[8].congestion_state)
	cong9.append(fac.stretches[0].cells[9].congestion_state)
	cong10.append(fac.stretches[0].cells[10].congestion_state)
	cong11.append(fac.stretches[0].cells[11].congestion_state)

	d.append(fac.stretches[0].delta_big[k-1])

print("\nEnd")
	
#print("Len rho: " + str(len(fac.stretches[0].cells[0].rho))) 

#############################
# Plot management section:  #
#############################
x_time = np.linspace(1, 24, 8640)
plt.figure(0)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('flow')
plt.plot(x_time, phi)

plt.figure(1)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r1')
plt.plot(x_time, r1)

plt.figure(2)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r2')
plt.plot(x_time, r2)

plt.figure(3)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r3')
plt.plot(x_time, r3)

plt.figure(4)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r4')
plt.plot(x_time, r4)

plt.figure(5)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r5')
plt.plot(x_time, r5)

plt.figure(6)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r6')
plt.plot(x_time, r6)

plt.figure(7)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r7')
plt.plot(x_time, r7)

plt.figure(8)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r8')
plt.plot(x_time, r8)

plt.figure(9)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r9')
plt.plot(x_time, r9)

plt.figure(10)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r10')
plt.plot(x_time, r10)

plt.figure(11)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r11')
plt.plot(x_time, r11)

plt.figure(12)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('r12')
plt.plot(x_time, r12)

plt.figure(90)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('fipiu')
plt.plot(x_time, fipiu)

plt.figure(91)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('fimeno')
plt.plot(x_time, fimeno)

plt.figure(80)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('e')
plt.plot(x_time, e0)

plt.figure(81)
plt.grid(True)
plt.xlabel('time [h]')
plt.ylabel('l')
plt.plot(x_time, l0)

plt.show()