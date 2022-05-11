import cell as c
import supervisor as s


sup = s.Supervisor()
sup.createCell(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
sup.createCell(10, 11, 12, 13, 14, 15, 16, 17, 18, 19)

print(sup.n_cells)
print()

for cell in sup.cells:
	cell.toString()
	print()