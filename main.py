import cell as c
import supervisor as s
import factory as f
import xlrd 

#loc = ("C:/A_Tesi/Python/CTM-s/CTM_data.xls")
loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")
wb = xlrd.open_workbook(loc)
sh = wb.sheet_by_index(0)

sh.cell_value(0,0)

fac = f.Factory()

fac.createStretch(10)

for i in range(1, sh.nrows):
	fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), 0, sh.cell_value(i,4), 0, 0, sh.cell_value(i,5), 0, 1)	

fac.addStationToStretch(0, 231, 1, 3, 890, 0.05, 0) #Note: q_max was statically assigned to the Qmax(4)/10 (the cell where the station merges back)

#for cell in fac.stretches[0].cells:
#	cell.toString()

#fac.stretches[0].stations[0].toString()
#fac.stretches[0].stations[0].services[0].toString()

k=0
while k<=1000:
	fac.stretches[0].update()
	print (k)
	k = k + 1
