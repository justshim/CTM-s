import csv

import numpy as np
from model import main_gdn
import pygad


if __name__ == '__main__':
    N_STRETCHES = 3000
    inputs = []
    phi = []

    with open('../data/phi_1_24h_realsmooth.csv', "r", encoding='utf8') as f_in:
        csv_reader = csv.reader(f_in)
        for row in csv_reader:
            phi.append(float(row[0]))

    for index in range(N_STRETCHES):
        print("\n")
        print("************************************")
        print("Stretch: " + str(index))
        print("************************************")
        print("\n")
        with open('../data/input1.csv', 'a', encoding='UTF8', newline='') as f_input:
            N_CELLS = 15
            stretch = np.zeros([N_CELLS, 6])
            for row in range(N_CELLS):
                stretch[row, 0] = row + 1
                stretch[row, 1] = np.random.randint(300, 1001)/1000  # L
                stretch[row, 2] = np.random.randint(80, 111)  # v
                stretch[row, 3] = np.random.randint(10, 41)  # w
                stretch[row, 4] = np.random.randint(1500, 2500)  # q_max
                stretch[row, 5] = np.random.randint(70, 100)  # rho_max
            np.savetxt(f_input, stretch, fmt='%g', delimiter=';')

            gen = 0

            def new_gen(ga):
                global gen

                gen = ga.generations_completed
                print("\n")
                print("************************************")
                print("gen: " + str(gen))
                print("************************************")
                print("\n")

            duration = 8640  # k=24h=8640 , k=1h=360, k=3h=1080
            rsmax = 1500

            station = [1, 1, 1]
            gene_type = [int, int, float]
            gene_space = [{'low': 1, 'high': N_CELLS-2}, {'low': 1, 'high': 720}, {'low': 0, 'high': 0.2}]

            parallel_processing = 16

            function_inputs = [station[0], station[1], station[2]]

            num_generations = 20
            num_parents_mating = 4

            sol_per_pop = 16
            num_genes = len(function_inputs)

            init_range_low = 0
            init_range_high = 720

            parent_selection_type = "sss"
            keep_parents = -1

            crossover_type = "two_points"
            crossover_probability = 0.5

            mutation_type = "random"
            mutation_num_genes = 1
            mutation_by_replacement = True
        #    initial_population = [[1, 2, 300, 0.1], [3, 6, 120, 0.05], [7, 9, 720, 0.2], [11, 12, 50, 0.02]]
            initial_population = None

            stop_criteria = "saturate_7"
            on_generation = new_gen

            with open('../data/ga_gdn1.csv', 'a', encoding='UTF8', newline='') as f:
                writer = csv.writer(f, delimiter=';')
                def fitness_func(solution, solution_idx):
                    output = [0, 0]
                    if solution[0] < solution[1]:
                        output = main_gdn.ga(stretch, phi, duration, rsmax, solution)
                        fit_int = 1000000 / output[0]
                        fit_pi = 20000 * output[1]
                        fitness = fit_int + fit_pi
                        print("Solution: [i: " + str(solution[0]) + ", j: " + str(solution[0]+2) +
                              ", delta: " + str(solution[1]) + ", beta: " + str(solution[2])+ "]\n" +
                              "Fitness Integral: " + str(fit_int)+
                              "\nFitness PI: " + str(fit_pi)+
                              "\nFitness: " + str(fitness))

                    else:
                        fitness = 0
                        print("Invalid solution: i>j")
                    return fitness

                fitness_function = fitness_func

                ga_instance = pygad.GA(num_generations=num_generations,
                                       num_parents_mating=num_parents_mating,
                                       fitness_func=fitness_function,
                                       sol_per_pop=sol_per_pop,
                                       num_genes=num_genes,
                                       gene_type=gene_type,
                                       gene_space=gene_space,
                                       on_generation=on_generation,
                                       stop_criteria=stop_criteria,
                                       init_range_low=init_range_low,
                                       init_range_high=init_range_high,
                                       mutation_by_replacement=mutation_by_replacement,
                                       initial_population=initial_population,
                                       crossover_probability=crossover_probability,
                                       parent_selection_type=parent_selection_type,
                                       keep_parents=keep_parents,
                                       crossover_type=crossover_type,
                                       mutation_type=mutation_type,
                                       mutation_num_genes=mutation_num_genes,
                                     #  mutation_percent_genes=mutation_percent_genes,
                                       parallel_processing=parallel_processing)

                ga_instance.run()

                solution, solution_fitness, solution_idx = ga_instance.best_solution()
                print("Parameters of the best solution : {solution}".format(solution=solution))
                print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
                writer.writerow([solution[0], solution[1], solution[2]])
                #ga_instance.plot_fitness(title="Fitness across generations")
