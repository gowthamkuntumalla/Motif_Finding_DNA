import os
from math import log
import time
import motif_finder
import matplotlib.pyplot as plt


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
        row_diff = sum((motif[i][j] * log(motif[i][j]/predictedmotif[i][j])  for j in range(len(motif[i]))))
        rel_ent += row_diff
    return rel_ent

def num_overlap_pos():
    """
    Number of overlapping positions between “sites.txt” and “predictedsites.txt”
    """
    pass


def num_overlap_sites(data_i = 0):
    """
    Number of overlapping sites (two sites overlap if at least ML/2 of their positions are common) between “sites.txt” and “predictedsites.txt
    ”"""
    sites = import_motif('results/dataset' + str(data_i) + '/sites.txt')
    predictedsites = import_motif('results/dataset' + str(data_i) + '/predictedsites.txt')


def runtime():
    """CPU runtime"""
    start = time.process_time() # CPU time
    ## code
    for i in range(10000):
        pass
    ##code
    end = time.process_time()
    print('CPU runtime: {0:.6f} sec'.format(end - start))

## Testing
KL_divergence()
