import cell as c
import supervisor as s
import factory as f
import xlrd 
import matplotlib.pyplot as plt
import numpy as np

#######################
# File read section:  #
#######################

## read file CTM_data from xls file

#loc = ("C:/A_Tesi/Python/CTM-s/CTM_data.xls")
loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")

wb = xlrd.open_workbook(loc)
sh = wb.sheet_by_index(0)
sh.cell_value(0,0)

## read phi first cell from xls file

#loc_phi = ("C:/A_Tesi/Python/CTM-s/phi_1.xls")
loc_phi = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/phi_1.xls")

wb_phi = xlrd.open_workbook(loc_phi)

sh_phi = wb_phi.sheet_by_index(0) 		#sheet 0 is a "realistic" input (24h)
#sh_phi = wb_phi.sheet_by_index(1)		#sheet 1 is a "synthetic" input with 2 equal peaks (24h)
#sh_phi = wb_phi.sheet_by_index(2)		#sheet 2 is a "synthetic" input with 1 peak (24h)
#sh_phi = wb_phi.sheet_by_index(3)		#sheet 3 is a "synthetic" input with 1 peak (3h)

sh_phi.cell_value(0,0)

phi_zero=[]
for i in range(0, sh_phi.nrows):
	phi_zero.append(sh_phi.cell_value(i,0))

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
	           #ID stretch   length,             v,                   w,              ,   q_max,           s,  r, rho_max,         beta, p
	fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), sh.cell_value(i,4), 0, 0, sh.cell_value(i,5), 0, 0.95)	
			#ID station, r_s_max, i, j, delta, beta_s, p

## create the stations via the factory
fac.addStationToStretch(0, 500, 3, 6, 5, 0.01, 1) #Note: r_s_max was statically assigned to the Qmax(4)/10 (the cell where the station merges back)

#for cell in fac.stretches[0].cells:
	#cell.toString()

#fac.stretches[0].stations[0].toString()

## support variables to save various parameters during execution, and possibly plot them
l0 = []

r0 = []
r1 = []
r2 = []
r3 = []
r4 = []
r5 = []
r6 = []
r7 = []
r8 = []

cong0 = []
cong1 = []
cong2 = []
cong3 = []
cong4 = []
cong5 = []
cong6 = []
cong7 = []
cong8 = []


##################################
# Exectution of the simulation:  #
##################################

k=0
while k<8640: 	# k=24h=8640 , k=1h=360, k=3h=1080
	print("Time instant: " + str(k) + "\n")
	
	fac.stretches[0].update(k)
	k = k + 1
	
	# save the various parameters in the previously created variables
	l0.append(fac.stretches[0].stations[0].l[k])
	
	r0.append(fac.stretches[0].cells[0].rho[k])
	r1.append(fac.stretches[0].cells[1].rho[k])
	r2.append(fac.stretches[0].cells[2].rho[k])
	r3.append(fac.stretches[0].cells[3].rho[k])
	r4.append(fac.stretches[0].cells[4].rho[k])
	r5.append(fac.stretches[0].cells[5].rho[k])
	r6.append(fac.stretches[0].cells[6].rho[k])
	r7.append(fac.stretches[0].cells[7].rho[k])
	r8.append(fac.stretches[0].cells[8].rho[k])

	cong0.append(fac.stretches[0].cells[0].congestionState[k])
	cong1.append(fac.stretches[0].cells[1].congestionState[k])
	cong2.append(fac.stretches[0].cells[2].congestionState[k])
	cong3.append(fac.stretches[0].cells[3].congestionState[k])
	cong4.append(fac.stretches[0].cells[4].congestionState[k])
	cong5.append(fac.stretches[0].cells[5].congestionState[k])
	cong6.append(fac.stretches[0].cells[6].congestionState[k])
	cong7.append(fac.stretches[0].cells[7].congestionState[k])
	cong8.append(fac.stretches[0].cells[8].congestionState[k])
	
#print("Len rho: " + str(len(fac.stretches[0].cells[0].rho))) 

#############################
# Plot management section:  #
#############################

x = np.linspace(0, 3, 8640)

plt.figure(0)
plt.grid(True)
plt.xlabel('hours')
plt.ylabel('rho')
plt.plot(x, l0)

# plt.figure(99)
# plt.grid(True)
# plt.plot(x, phi_zero)
#plt.scatter(np.linspace(0, 8640, 8640), r0, 4, marker="x")


# plt.figure(1)
# plt.grid(True)
# plt.plot(r1)

# plt.figure(2)
# plt.grid(True)
# plt.plot(r2)

# plt.figure(3)
# plt.grid(True)
# plt.plot(r3)

# plt.figure(4)
# plt.grid(True)
# plt.plot(r4)

# plt.figure(5)
# plt.grid(True)
# plt.plot(r5)

# plt.figure(6)
# plt.grid(True)
# plt.plot(r6)

# plt.figure(7)
# plt.grid(True)
# plt.plot(r7)

# plt.figure(8)
# plt.grid(True)
# plt.plot(r8)

# x = np.linspace(0, 8640, 8640)
# y = r0
# z = np.polyfit(x, y, 10)
# p = np.poly1d(z)

#plt.plot(x,p(x),"r--")

plt.show()
