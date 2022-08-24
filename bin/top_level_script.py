import time

import numpy

from model import main_ga
import pygad

# from model.optimization import data_creation as opt
#
# if __name__ == '__main__':
#     loc = ("H:/Il mio Drive/Tesi magistrale/CTMs-identification/fnc/extracted_data/CTM_param_out_nice.xls")
#     # loc = ("C:/A_Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls")
#     # loc = "C:/Users/adria/Documents/Uni/LM II anno/Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out.xls"
#     path_file_output = "../data/opti_data.csv"
#     duration = 8640  # k=24h=8640 , k=1h=360, k=3h=1080
#     t = time.time()
#     o = opt.data_creation(loc, path_file_output, duration)
#     o.initVars()
#     o.firstIteration()
#     # initDelta, delta, stepDelta, initBeta, beta, stepBeta, initPriority (MS), priority (MS), stepPriority (MS)
#     o.initSimulation(60, 721, 60, 1, 21, 1, 95, 96, 1)
#     elapsed = time.time() - t
#     print("Elapsed time: " + str(elapsed))


if __name__ == '__main__':
    t = time.time()

    # path = ("H:/Il mio Drive/Tesi magistrale/CTMs-identification/fnc/extracted_data/CTM_param_out_nice.xls")
    # path = ("C:/A_Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out_nice.xls")
    path = "C:/Users/adria/Documents/Uni/LM II anno/Tesi/CTMs-identification/fnc/extracted_data/CTM_param_out_nice.xls"

    path_file_output = "../data/ga_result.csv"

    duration = 8640  # k=24h=8640 , k=1h=360, k=3h=1080
    # onramps = [[1000, 3000, 4, 0.01]]
    # offramps = [[2, 0.5]]
    onramps = []
    offramps = []

    station = [1, 1, 1, 1, 1, 1]

    gene_type = [int, int, int, int, float, float]
    gene_space = [{'low': 500, 'high': 500}, {'low': 1, 'high': 12}, {'low': 2, 'high': 13},
                  {'low': 60, 'high': 720}, {'low': 0, 'high': 0.2}, {'low': 0, 'high': 0.1}]

    parallel_processing = 4

    function_inputs = [station[0], station[1], station[2], station[3], station[4], station[5]]
    desired_output = 0

    num_generations = 50
    num_parents_mating = 4

    sol_per_pop = 8
    num_genes = len(function_inputs)

    init_range_low = -2
    init_range_high = 5

    parent_selection_type = "sss"
    keep_parents = 1

    crossover_type = "single_point"

    mutation_type = "random"
    mutation_percent_genes = 50
    initial_population = [[500, 1, 2, 500, 0.1, 0.05], [500, 2, 3, 500, 0.2, 0.05], [500, 1, 2, 500, 0.1, 0.05], [500, 2, 3, 500, 0.2, 0.05]]


    def fitness_func(solution, solution_idx):
        output = [0, 0]
        if solution[1] < solution[2]:
            output = main_ga.ga(path, duration, solution, onramps, offramps)
            fitness = 10000 / numpy.abs(output[0] - desired_output)
            print("Solution: [rs_max: " + str(solution[0]) + ", i: " + str(solution[1]) + ", j: " + str(solution[2]) +
                  ", delta: " + str(solution[3]) + ", beta: " + str(solution[4]) + ", p: " + str(solution[5]) + "]\n" +
                  "Fitness: " + str(fitness))
        else:
            fitness = -99999999
            print("Invalid solution mi vergogno sono un incapace cane randagio")
        return fitness


    fitness_function = fitness_func

    ga_instance = pygad.GA(num_generations=num_generations,
                           num_parents_mating=num_parents_mating,
                           fitness_func=fitness_function,
                           sol_per_pop=sol_per_pop,
                           num_genes=num_genes,
                           gene_type=gene_type,
                           gene_space=gene_space,
                           init_range_low=init_range_low,
                           init_range_high=init_range_high,
                           initial_population=initial_population,
                           parent_selection_type=parent_selection_type,
                           keep_parents=keep_parents,
                           crossover_type=crossover_type,
                           mutation_type=mutation_type,
                           mutation_percent_genes=mutation_percent_genes,
                           parallel_processing=parallel_processing)

    ga_instance.run()

    solution, solution_fitness, solution_idx = ga_instance.best_solution()
    print("Parameters of the best solution : {solution}".format(solution=solution))
    print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))

    prediction = numpy.sum(numpy.array(function_inputs) * solution)
    print("Predicted output based on the best solution : {prediction}".format(prediction=prediction))

    elapsed = time.time() - t
    print("Elapsed time: " + str(elapsed))
