import cell as c
import supervisor as s
import factory as f
import xlrd 
import matplotlib.pyplot as plt
import numpy as np

loc = ("C:/A_Tesi/Python/CTM-s/CTM_data.xls")
loc_phi = ("C:/A_Tesi/Python/CTM-s/phi_1.xls")
#loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")
wb = xlrd.open_workbook(loc)
sh = wb.sheet_by_index(0)

wb_phi = xlrd.open_workbook(loc_phi)
sh_phi = wb_phi.sheet_by_index(0)

sh.cell_value(0,0)
sh_phi.cell_value(0,0)

phi_zero=[]
for i in range(0, sh_phi.nrows):
	phi_zero.append(sh_phi.cell_value(i,0))

fac = f.Factory()

### T[h]
fac.createStretch(10/3600, 5000, phi_zero)

for i in range(1, sh.nrows):
	           #ID stretch   length,             v,                   w,                q,   q_max,               s,    r, rho_max,          beta, p
	fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), 1500, sh.cell_value(i,4), 2500, 0, sh.cell_value(i,5), 0, 0.99)	

fac.addStationToStretch(0, 231, 1, 3, 890, 0.05, 1) #Note: q_max was statically assigned to the Qmax(4)/10 (the cell where the station merges back)

#for cell in fac.stretches[0].cells:
	#cell.toString()

#fac.stretches[0].stations[0].toString()

k=0

r0 = []
r1 = []
r2 = []
r3 = []
r4 = []
r5 = []
r6 = []
r7 = []
r8 = []

## k=8640=24h
while k<=8639:
	print()
	#print("Time instant: " + str(k))
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
	
plt.figure(0)
plt.plot(r0)

plt.figure(1)
plt.plot(r1)

plt.figure(2)
plt.plot(r2)

plt.figure(3)
plt.plot(r3)

plt.figure(4)
plt.plot(r4)

plt.figure(5)
plt.plot(r5)

plt.figure(6)
plt.plot(r6)

plt.figure(7)
plt.plot(r7)

plt.figure(8)
plt.plot(r8)

plt.show()



# ATTENZIONE Ss SEMPRE 0