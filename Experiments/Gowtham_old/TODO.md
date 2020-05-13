1. Finish the normalization functions in motif_finder.py (pg 78, 79 of bailey elkan 1995)

2. verify all the functinality in evaluate.py and motif_finder.py
predicted motif.txt seems to be wrong. check that carefully.


Imporvements (try definitely):
1. Vary the n_iter from [10, 100, 1000, 10000] - 1000 takes 6 hours per run (all 70 datasets)

2. Change the starting value of priors (latent variables) as per Bailey Elkan 1994.pdf Implementation section.

3. Try changing site prediction criterion from max(row) to something else. 
