import cell as c
import supervisor as s
import factory as f
import xlrd 

loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")
wb = xlrd.open_workbook(loc)
sh = wb.sheet_by_index(0)

sh.cell_value(0,0)

fac = f.Factory()

fac.createStretch(10)

for i in range(1, sh.nrows):
	fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), 0, sh.cell_value(i,4), 0, 0, sh.cell_value(i,5), 0, 0)	

fac.addStationToStretch(0, 231, 1, 3) #Note: q_max was statically assigned to the Qmax(4)/10 (the cell where the station merges back)
fac.addServiceToStation(0, 0, 890, 0.05)

for cell in fac.stretches[0].cells:
	cell.toString()

fac.stretches[0].stations[0].toString()
fac.stretches[0].stations[0].services[0].toString()