import os
import sys
from math import log
import numpy as np
from scipy.special import rel_entr
import time
from motif_finder import MotifFinder
import matplotlib.pyplot as plt
from shutil import copytree,copyfile # for testing

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


"""
Testing

Iterate over parameter. Average for i to i+9 for each parameter set

Generating parameters set (length = 7) for iteration. 
See benchmarks.py for order. 
Each parameter_set has 10 datasets in that order

"""

icpc = [1.0, 1.5, 2.0] # info_cont_per_col
ml = [6, 7, 8] # motif_length
sl = 500 # sequence_length
sc = [5, 10, 20] #sequence_count
nucleotides = ['A', 'G', 'C', 'T']
default_ml = 7
default_icpc = 1.5
default_sc = 10

# Populating parameters set
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

try:
	os.mkdir('results')

except FileExistsError:
	print('folder results: already exists')

for i,pset in enumerate(params_set):
	kl_div = []
	overlap_pos = [] 
	overlap_sites = [] 
	runtimes = []

	for j in range(10): 
	# number of repetitions to get a statistical average
		filenum = i*10+j
		
		try:
			src_folder = 'benchmarks/dataset' + str(filenum)
			dest_folder = 'results/dataset' + str(filenum)
			copytree(src_folder,dest_folder)

			## ***** copy for testing. change it after implementing motif_finder.py *****
			copyfile(dest_folder + '/sites.txt', dest_folder + '/predictedsites.txt')
			copyfile(dest_folder + '/motif.txt', dest_folder + '/predictedmotif.txt')
		
		except:
			# print('folder already exists',sys.exc_info()[0])
			pass

		"""

		Calling main algorithm

		"""
		sol = MotifFinder()
		sol.set_motif_length(dest_folder + '/motiflength.txt')
		sol.import_sequences(dest_folder + '/sequences.fa')

		start_time = time.process_time() # CPU time: run start

		sol.predict_motif() # creates predictedmotif.txt
		sol.predict_sites() # creates predictedsites.txt

		end_time = time.process_time() # CPU time: run end 

		"""
		
		Calling evaluator functions

		"""

		kl_div.append(KL_divergence(filenum))
		overlap_pos.append(num_overlap_pos(filenum))
		overlap_sites.append(num_overlap_sites(filenum))
		runtimes.append(round(end_time - start_time,8))

	result_kl_div.append(np.mean(kl_div))
	result_overlap_pos.append(np.mean(overlap_pos)) 
	result_overlap_sites.append(np.mean(overlap_sites))
	result_runtime.append(np.mean(runtimes))


print("KL_Divergence:",result_kl_div) 
print("Overlapping Postions",result_overlap_pos) 
print("Overlapping Sites:",result_overlap_sites)
print("Runtimes",result_runtime)

try:
	os.mkdir('performance_plots')

except FileExistsError:
	print('folder performace_plots: already exists')

"""	
Plotting Results

"""

def plotter(X,Y, y_lab = None):
	plt.plot(np.arange(len(params_set)), Y)
	plt.title(y_lab + ' vs Parameter Set')
	plt.xlabel('Parameter Set')
	plt.ylabel(y_lab)
	plt.ticklabel_format(axis = 'y',style = 'sci')

	for i, txt in enumerate(X):
		plt.annotate(tuple(txt), (i, Y[i]))

	plt.savefig('performance_plots/'+ y_lab +'.png', dpi = 150)
	plt.close()

plotter(params_set, result_kl_div, y_lab = 'Relative Entropy')
plotter(params_set, result_overlap_pos, y_lab = 'Overlapping Positions')
plotter(params_set, result_overlap_sites, y_lab = 'Overlapping Sites')
plotter(params_set, result_runtime, y_lab = 'Runtime')

