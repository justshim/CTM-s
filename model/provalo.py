import itertools
import numpy as np
import pandas as pd

if __name__ == '__main__':
    generation = False
    #print("GPU test: " + str(tf.test.is_gpu_available()) + "\n")
    path_file_output = "../data/saaaaaaa.csv"
    if(generation == True):
        print("Start: ")
        aa = 16
        bb = aa+1
        a = range(aa, aa+1, 1)      #i
        b = range(bb, bb+4, 1)      #j
        c = range(90, 721, 90)  #delta
        d = range(1, 21, 2)     #beta
        e = range(1, 21, 2)     #priority
        f = range(100, 1001, 200)      #L
        g = range(80, 121, 10)      #v
        h = range(10, 51, 10)      #w
        i = range(1500, 2501, 200)  #qmax
        l = range(50, 101, 10)     #rhomax

        #paramlist = list(itertools.product(a, b, c, d, e, f, g, h, i, l, ff, gg, hh, ii, ll))  # all possible combinations of tuples
        paramlist = list(itertools.product(a, b, c, d, e, f, g, h, i, l))  # all possible combinations of tuples
        print("Total cases to be evaluated: " + str(len(paramlist)))

        arr = np.array(paramlist)
        with open(path_file_output, 'a') as f:
            np.savetxt(f, arr, delimiter=";", fmt='%g')
    if(generation == False):


        noise = np.random.normal(100, 100, 100)
        print(noise)
        #'L_i','v_i','w_i','qmax_i','rhomax_i'
        fields = ['L_i']
        #a = pd.read_csv("C:/Users/dspal/Desktop/super_mega_enorme_dataset.csv", usecols=fields, sep=';')

        # See content in 'star_name'
        print("ghe semo: ")





