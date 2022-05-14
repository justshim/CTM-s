import cell as c
import supervisor as s
import factory as f
import xlrd 
import matplotlib.pyplot as plt
import numpy as np

#loc = ("C:/A_Tesi/Python/CTM-s/CTM_data.xls")
loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")
wb = xlrd.open_workbook(loc)
sh = wb.sheet_by_index(0)

sh.cell_value(0,0)

fac = f.Factory()

fac.createStretch(10, 0, 0)

for i in range(1, sh.nrows):
	           #ID stretch   length,             v,                   w,                q,   q_max,            s, r, rho_max,          beta, p
	fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), 1, sh.cell_value(i,4), 1, 1, sh.cell_value(i,5), 1, 0.8)	

fac.addStationToStretch(0, 231, 1, 3, 890, 0.05, 1) #Note: q_max was statically assigned to the Qmax(4)/10 (the cell where the station merges back)

#for cell in fac.stretches[0].cells:
	#cell.toString()

#fac.stretches[0].stations[0].toString()

k=0
r = []
while k<=1000:
	fac.stretches[0].update()
	#print (k)
	#print(fac.stretches[0].cells[-1].rho)
	k = k + 1
	#r.append(fac.stretches[0].cells[-2].phi)

#plt.plot(r)
#plt.show()


