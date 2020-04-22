import os
from math import log
import numpy as np
from scipy.special import rel_entr
import time
import motif_finder
import matplotlib.pyplot as plt
from shutil import copytree,copyfile # for testing

try:
    os.mkdir('results')
except FileExistsError:
    pass
    #print('folder results: already exists')

def import_motif(file):
    """
    Used by function KL_divergence
        :file = complete file location + motif file name
    """
    with open(file,'r') as f:
        lines = f.readlines()
    motif = []
    for line in lines[1:-1]:
        motif.append(list(map(float,line.strip().split())))
    return motif
    

def KL_divergence(data_i = 0):
    """
    Relative Entropy b/w
    
    P: True distribution = “motif.txt”
        and
    Q: Model Predicted/Approximate distribution = “predictedmotif.txt”
    
    Output:
        : D_KL(P||Q)
    """
    motif = import_motif('results/dataset' + str(data_i) + '/motif.txt') # list of lists
    predictedmotif = import_motif('results/dataset' + str(data_i) + '/predictedmotif.txt')
    
    rel_ent = 0
    
    for i in range(len(motif[1:-1])):
        # compare each row (ACGT) against each other in two matrices
        # row_diff = sum((motif[i][j] * log(motif[i][j]/predictedmotif[i][j])  for j in range(len(motif[i]))))
        row_diff = rel_entr(motif[i],predictedmotif[i])
        rel_ent += row_diff
    return rel_ent

def num_overlap_pos(data_i = 0):
    """
    Number of overlapping positions between “sites.txt” and “predictedsites.txt”
    """
    result = 0
    return result


def num_overlap_sites(data_i = 0):
    """
    Number of overlapping sites (two sites overlap if at least ML/2 of their positions are common) between “sites.txt” and “predictedsites.txt
    ”"""
    sites = import_motif('results/dataset' + str(data_i) + '/sites.txt')
    predictedsites = import_motif('results/dataset' + str(data_i) + '/predictedsites.txt')

    result = 0
    return result


def runtime(data_i = 0):
    """CPU runtime"""
    start = time.process_time() # CPU time
    ## code
    for i in range(10000):
        pass
    ##code
    end = time.process_time()
    # print('CPU runtime: {0:.6f} sec'.format(end - start))

    return round(end - start,6)


"""
Testing

Iterate over parameter. Average for i to i+9 for each parameter set
"""

icpc = [1.0, 1.5, 2.0] # info_cont_per_col
ml = [6, 7, 8] # motif_length
sl = 500 # sequence_length
sc = [5, 10, 20] #sequence_count
nucleotides = ['A', 'G', 'C', 'T']
default_ml = 7
default_icpc = 1.5
default_sc = 10


""" Generating parameters set (length = 7) for iteration. See benchmarks.py for order. Each parameter_set has 10 datasets in that order """  
params_set = [] 

for i in icpc: # count = 3
	params_set.append([i,default_ml,sl,default_sc]) 

for i in [6,8]: # count = 2
	params_set.append([default_icpc,i,sl,default_sc]) 

for i in [5,20]: # count = 2
	params_set.append([default_icpc,default_ml,sl,i]) 


## Result arrays. Length = 7
result_kl_div = [] 
result_overlap_pos = [] 
result_overlap_sites = [] 
result_runtime = []

for i,pset in enumerate(params_set):
	kl_div = []
	overlap_pos = [] 
	overlap_sites = [] 
	runtimes = []
	for j in range(10):
		filenum = i*10+j
		
		try:
			copytree('benchmarks/dataset'+ str(data_i),'results/dataset' + str(data_i))
			## ***** copy for testing. change it after implementing motif_finder.txt *****
			copyfile('results/dataset' + str(data_i) + '/sites.txt',  'results/dataset' + str(data_i) + '/predictedsites.txt')
			copyfile('results/dataset' + str(data_i) + '/motif.txt',  'results/dataset' + str(data_i) + '/predictedmotif.txt')
		except:
			pass

		kl_div.append(KL_divergence(filenum))
		overlap_pos.append(num_overlap_pos(filenum))
		overlap_sites.append(num_overlap_sites(filenum))
		runtimes.append(runtime(filenum))

	result_kl_div.append(np.mean(kl_div))
	result_overlap_pos.append(np.mean(overlap_pos)) 
	result_overlap_sites.append(np.mean(overlap_sites))
	result_runtime.append(np.mean(runtimes))


print(result_kl_div) 
print(result_overlap_pos) 
print(result_overlap_sites)
print(result_runtime)

try:
	os.mkdir('performance_plots')

except FileExistsError:
    pass
    #print('folder results: already exists')


def plotter(X,Y, y_lab = None):
	plt.plot(np.arange(len(params_set)), Y)
	plt.title(y_lab + ' vs Parameter Set')
	plt.xlabel('Parameter Set')
	plt.ylabel(y_lab)
	plt.ticklabel_format(axis = 'y',style = 'sci')
	for i, txt in enumerate(X):
	    plt.annotate(tuple(txt), (i, Y[i]))
	    # plt.annotate(i, (i, Y[i]))

	plt.savefig('performance_plots/'+ y_lab +'.png', dpi = 150)
	plt.close()

plotter(params_set, result_kl_div, y_lab = 'Relative Entropy')
plotter(params_set, result_overlap_pos, y_lab = 'Overlapping Positions')
plotter(params_set, result_overlap_sites, y_lab = 'Overlapping Sites')
plotter(params_set, result_runtime, y_lab = 'Runtime')

