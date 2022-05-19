import cell as c
import supervisor as s
import factory as f
import xlrd 
import matplotlib.pyplot as plt
import numpy as np

# read file CTM_data
loc = ("C:/A_Tesi/Python/CTM-s/CTM_data.xls")
#loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")
wb = xlrd.open_workbook(loc)
sh = wb.sheet_by_index(0)
sh.cell_value(0,0)

#read phi first cell 
loc_phi = ("C:/A_Tesi/Python/CTM-s/phi_1.xls")
#loc_phi = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/phi_1.xls")

wb_phi = xlrd.open_workbook(loc_phi)

#sh_phi = wb_phi.sheet_by_index(0) 		#sheet 0 is a "realistic" input
#sh_phi = wb_phi.sheet_by_index(1)		#sheet 1 is a "synthetic" input with 2 equal peaks
#sh_phi = wb_phi.sheet_by_index(2)		#sheet 2 is a "synthetic" input with 1 peak (24h)
sh_phi = wb_phi.sheet_by_index(3)		#sheet 3 is a "synthetic" input with 1 peak (3h)

sh_phi.cell_value(0,0)

phi_zero=[]
for i in range(0, sh_phi.nrows):
	phi_zero.append(sh_phi.cell_value(i,0))

fac = f.Factory()

	# timeLength [h],   lastPhi,  phi_zero
fac.createStretch(10/3600, 5000, phi_zero) 

for i in range(1, sh.nrows):
	           #ID stretch   length,             v,                   w,              ,   q_max,               s,  r, rho_max,         beta, p
	fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), sh.cell_value(i,4), 2500, 0, sh.cell_value(i,5), 0, 0.95)	
			#ID station, r_s_max, i, j, delta, beta_s, p
fac.addStationToStretch(0, 231, 3, 6, 890, 0.05, 1) #Note: q_max was statically assigned to the Qmax(4)/10 (the cell where the station merges back)

#for cell in fac.stretches[0].cells:
	#cell.toString()

#fac.stretches[0].stations[0].toString()

r0 = []
r1 = []
r2 = []
r3 = []
r4 = []
r5 = []
r6 = []
r7 = []
r8 = []

k=0
while k<1080: 						# k=24h=8640 , k=1h=360, k=3h=1080
	#print("Time instant: " + str(k) + "\n")
	
	fac.stretches[0].update(k)
	k = k + 1
	r0.append(fac.stretches[0].cells[0].rho[k])
	r1.append(fac.stretches[0].cells[1].rho[k])
	r2.append(fac.stretches[0].cells[2].rho[k])
	r3.append(fac.stretches[0].cells[3].rho[k])
	r4.append(fac.stretches[0].cells[4].rho[k])
	r5.append(fac.stretches[0].cells[5].rho[k])
	r6.append(fac.stretches[0].cells[6].rho[k])
	r7.append(fac.stretches[0].cells[7].rho[k])
	r8.append(fac.stretches[0].cells[8].rho[k])
	
#print("Len rho: " + str(len(fac.stretches[0].cells[0].rho))) 

x = np.linspace(0, 3, 1080)

plt.figure(0)
plt.grid(True)
plt.plot(x, r7)

plt.figure(99)
plt.grid(True)
plt.plot(x, phi_zero)
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
