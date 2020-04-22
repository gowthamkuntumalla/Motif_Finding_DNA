import os
import time
import motif_finder
import matplotlib.pyplot as plt


try:
    os.mkdir('results')
except FileExistsError:
    print('folder results: already exists')

def KL_divergence():
    """
    Relative Entropy between “motif.txt” and “predictedmotif.txt”
    """
    pass


def num_overlap_pos():
    """
    Number of overlapping positions between “sites.txt” and “predictedsites.txt”
    """
    pass


def num_overlap_sites():
    """
    Number of overlapping sites (two sites overlap if at least ML/2 of their positions are common) between “sites.txt” and “predictedsites.txt
    ”"""
    pass


def runtime():
    """CPU runtime"""
    start = time.process_time() # CPU time
    
    ## code
    for i in range(10000):
        pass
    ##code
    
    end = time.process_time()
    
    print('CPU runtime: {0:.6f} sec'.format(end - start))

runtime()
