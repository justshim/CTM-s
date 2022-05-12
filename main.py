import cell as c
import supervisor as s


stretch2 = s.Stretch(1)
stretch2.createCell(20, 21, 22, 23, 24, 25, 26, 27, 28, 29)
stretch2.createCell(30, 31, 32, 33, 34, 35, 36, 37, 38, 39)
stretch2.createCell(40, 41, 42, 43, 44, 45, 46, 47, 48, 49)

stretch2.createStation(1, 2, 3)
stretch2.createStation(11, 22, 33)

stretch2.stations[0].createService(1,2)
stretch2.stations[0].createService(88,888)
stretch2.stations[1].createService(99,999)

# print("**********\nStretch 2\n**********")
#print(stretch2.n_cells)
#print()
stretch2.update()

for cell in stretch2.cells:
	cell.toString()
	print()

#print("====================")

#for station in stretch2.stations:
#	station.toString()
#	print()	

#print("====================")

#for station in stretch2.stations:
#	for i in range(len(station.services)):
#		station.services[i].toString()
#		print()

