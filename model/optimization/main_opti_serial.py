from model import cell as c
from model import supervisor as s
from model import factory as f
import xlrd
import xlwt
import matplotlib.pyplot as plt
import numpy as np

#######################
# File read section:  #
#######################
#path_file_output = 'C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/opti_data.xls'
path_file_output = 'opti_data.xls'


## read file CTM_data from xls file

#loc = ("C:/A_Tesi/Python/CTM-s/CTM_data.xls")
#loc = ("C:/A_Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
loc = ("H:/Il mio Drive/Tesi magistrale/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
#loc = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/CTM_data.xls")

wb = xlrd.open_workbook(loc)
sh = wb.sheet_by_name("Cells parameters")
sh.cell_value(0,0)

wbt = xlwt.Workbook()
ws = wbt.add_sheet('Sheet1')

## read phi first cell from xls file

#loc_phi = ("C:/A_Tesi/Python/CTM-s/phi_1.xls")
#loc_phi = ("C:/A_Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
#loc_phi = ("C:/Users/adria/Documents/Uni/LM II anno/Tesi/python/CTM-s/phi_1.xls")
loc_phi=loc
wb_phi = xlrd.open_workbook(loc_phi)

#sh_phi = wb_phi.sheet_by_index(0)         #sheet 0 is a "realistic" input (24h)
#sh_phi = wb_phi.sheet_by_index(1)        #sheet 1 is a "synthetic" input with 2 equal peaks (24h)
#sh_phi = wb_phi.sheet_by_index(2)        #sheet 2 is a "synthetic" input with 1 peak (24h)
#sh_phi = wb_phi.sheet_by_index(3)        #sheet 3 is a "synthetic" input with 1 peak (3h)
#sh_phi = wb_phi.sheet_by_index(4)        #sheet 4 is a flat input (24h)

sh_phi = wb_phi.sheet_by_name("First Demand Smooth")
sh_phi.cell_value(0,0)
phi_zero=[]
for i in range(0, sh_phi.nrows):
    phi_zero.append(sh_phi.cell_value(i,0))

sh_last_phi = wb.sheet_by_name("Last Demand Smooth")
sh_last_phi.cell_value(0,0)
last_phi=[]
for i in range(0, sh_last_phi.nrows):
    last_phi.append(sh_last_phi.cell_value(i,0))

cell_in = 0
row = 0
duration = 8640 # k=24h=8640 , k=1h=360, k=3h=1080

####################################
# First iteration without stations #
####################################
fac = f.Factory()

## create the stretch via the factory
    # timeLength [h],   lastPhi,  phi_zero
fac.createStretch(sh.cell_value(2,6), last_phi, phi_zero)

## create the cells via the factory
for i in range(1, sh.nrows):
               #ID stretch   length, v_free,                   w,              ,   q_max, rho_max,         p_ms
    fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), sh.cell_value(i,4), sh.cell_value(i,5), 0.95)

## create the on-ramps via the factory
            #ID stretch, d_r, r_r_max, j, p_r
#fac.addOnRampToStretch(0, 100, 500, 2, 0.05)

## create the off-ramps via the factory
            #ID_stretch, i, beta_r
#fac.addOffRampToStretch(0, 7, 0.05)

ttt = fac.stretches[0].computeTTT()

d0 = []

##################################
# Exectution of the simulation:  #
##################################

k=0
while k<duration:

    fac.stretches[0].update(k)
    k = k + 1

    d0.append(fac.stretches[0].delta_big[k-1])

delta_max0 = max(d0)
print("max(delta_0): " + str(delta_max0))
print()

ws.write(0, 0, "i")
ws.write(0, 1, "j")
ws.write(0, 2, "delta")
ws.write(0, 3, "beta")
ws.write(0, 4, "integral")
ws.write(0, 5, "max_delta")
ws.write(0, 6, "pi")



for cell_in in range(1, sh.nrows-1, 1):

    for cell_out in range(cell_in+1, sh.nrows, 1):

        for delta in range(60, 61, 60): #360k = 3600s = 60 min. 720k = 7200s = 120 min

            for beta in range(1, 2, 1):
                beta=beta/100
###################################################
                # Initialization of all components of the model: #
###################################################

                ## create a factory instance that manages the creation of objects
                fac = f.Factory()

                ## create the stretch via the factory
                    # timeLength [h],   lastPhi,  phi_zero
                fac.createStretch(sh.cell_value(2,6), last_phi, phi_zero)

                ## create the cells via the factory
                for i in range(1, sh.nrows):
                               #ID stretch   length, v_free,                   w,              ,   q_max, rho_max,         p_ms
                    fac.addCellToStretch(0, sh.cell_value(i,1), sh.cell_value(i,2), sh.cell_value(i,3), sh.cell_value(i,4), sh.cell_value(i,5), 0.95)

                ## create the stations via the factory
                            #ID stretch, r_s_max, i, j, delta, beta_s, p
                fac.addStationToStretch(0, 500, cell_in, cell_out, delta, beta, 0.05)

                ## create the on-ramps via the factory
                            #ID stretch, d_r, r_r_max, j, p_r
                #fac.addOnRampToStretch(0, 100, 500, 2, 0.05)

                ## create the off-ramps via the factory
                            #ID_stretch, i, beta_r
                #fac.addOffRampToStretch(0, 7, 0.05)

                #ttt = fac.stretches[0].computeTTT()


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

                cong0 = []
                cong1 = []
                cong2 = []
                cong3 = []
                cong4 = []
                cong5 = []
                cong6 = []
                cong7 = []
                cong8 = []

                d = []

                ##################################
                # Exectution of the simulation:  #
                ##################################

                k=0
                while k<duration:
                    #print("\nTime instant: " + str(k))

                    fac.stretches[0].update(k)
                    k = k + 1

                    # save the various parameters in the previously created variables
                    #l0.append(fac.stretches[0].stations[0].l[k])
                    #e0.append(fac.stretches[0].stations[0].e[k])

                    r0.append(fac.stretches[0].cells[0].rho[k])
                    r1.append(fac.stretches[0].cells[1].rho[k])
                    r2.append(fac.stretches[0].cells[2].rho[k])
                    r3.append(fac.stretches[0].cells[3].rho[k])
                    r4.append(fac.stretches[0].cells[4].rho[k])
                    r5.append(fac.stretches[0].cells[5].rho[k])
                    r6.append(fac.stretches[0].cells[6].rho[k])
                    r7.append(fac.stretches[0].cells[7].rho[k])
                    r8.append(fac.stretches[0].cells[8].rho[k])

                    cong0.append(fac.stretches[0].cells[0].congestion_state)
                    cong1.append(fac.stretches[0].cells[1].congestion_state)
                    cong2.append(fac.stretches[0].cells[2].congestion_state)
                    cong3.append(fac.stretches[0].cells[3].congestion_state)
                    cong4.append(fac.stretches[0].cells[4].congestion_state)
                    cong5.append(fac.stretches[0].cells[5].congestion_state)
                    cong6.append(fac.stretches[0].cells[6].congestion_state)
                    cong7.append(fac.stretches[0].cells[7].congestion_state)
                    cong8.append(fac.stretches[0].cells[8].congestion_state)

                    d.append(fac.stretches[0].delta_big[k-1])

                integ = 0
                for kkk in d:
                    integ = integ + kkk

                print("i: " + str(cell_in) + "; j: " + str(cell_out))
                print("delta stazione: " + str(delta))
                print("beta: " + str(beta))
                print("Integ: " + str(integ))
                delta_max = max(d)
                pi = (delta_max0-delta_max)/delta_max0
                print("Pi: " + str(pi))
                print()

                row = row + 1
                ws.write(row, 0, cell_in)
                ws.write(row, 1, cell_out)
                ws.write(row, 2, delta)
                ws.write(row, 3, beta)
                ws.write(row, 4, integ)
                ws.write(row, 5, delta_max)
                ws.write(row, 6, pi)

    #     plt.figure(row)
    #     plt.grid(True)
    #     plt.plot(cong8)

    # plt.show()

#wbt.save(path_file_output)
wbt.save(path_file_output)
#print("Len rho: " + str(len(fac.stretches[0].cells[0].rho)))

#############################
# Plot management section:  #
#############################

# plt.figure(0)
# plt.grid(True)
# plt.xlabel('k')
# plt.ylabel('r6')
# plt.plot(d)



