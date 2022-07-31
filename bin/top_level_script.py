#from model import main as model
import time
from model.optimization import data_creation as opt

if __name__ == '__main__':
	loc = ("H:/Il mio Drive/Tesi magistrale/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
	path_file_output = ("../data/opti_data.csv")
	duration = 8640 # k=24h=8640 , k=1h=360, k=3h=1080
	t = time.time()
	o = opt.data_creation(loc, path_file_output, duration)
	o.initVars()
	o.firstIteration()
	#initDelta, delta, stepDelta, initBeta, beta, stepBeta, initPriority, priority, stepPriority
	o.initSimulation(60, 61, 60, 7, 8, 1, 94, 95, 1)
	elapsed = time.time() - t
	print("Elapsed time: " + str(elapsed))