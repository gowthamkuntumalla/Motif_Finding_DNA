#!/usr/bin/env python3

""" Algorithm to predict motif (position weight matrix) and predict sites

	Called:
	By evaluate.py
	
	Input
	: In folder 'results/dataset{i}'
	:: motiflength.txt
	:: sequences.fa
	
	Output
	: In folder 'results/dataset{i}'
	:: predictedmotif.txt
	:: predictedsites.txt
"""


import os
import numpy as np
from functools import reduce 
import operator

class MotifFinder:
	""" Solution Class """
	
	def __init__(self):
		self.sequenceList = []
		self.motifLen = 0
		self.n_seq = 0 # number of sequences in the file sequences.fa
		self.seqLen = 0 # avg length of a sequence.
		self.alphabet = {'A':0, 'G':1, 'C':2, 'T':3} # positions of Nucelic Acid alphabet
		self.lenAlphabet = len(self.alphabet.keys())
	
	def set_motif_length(self,lengthfile):
		""" 
		number of DNA sequences in the file sequences.fa 
		"""
		with open(lengthfile,'r') as f:
			lines = f.readlines()
		self.motifLen = int(lines[0])
	
	
	def import_sequences(self,sequencesFile):
		""" 
		sequencesFile is the file name in results/dataset{i}
		"""
		with open(sequencesFile,'r') as f:
			lines = f.readlines()

		for i,line in enumerate(lines):
			if i%2 == 1:
				self.sequenceList.append(line.strip())

		self.n_seq = len(self.sequenceList)
		self.seqLen = int(np.mean([len(l) for l in self.sequenceList]))

	
	
	def _conditional_proba(self,seq,theta0,theta1):
		"""
		Private func.
		To be called in expectation step of iterations

		output:
			tuple(P(Xi|theta1),P(Xi|theta2))
		"""

		p_x_t1 = 1 # motif class
		p_x_t2 = 1 # background class
	
		for i,c in enumerate(seq):
			p_x_t1 *= theta1[i][self.alphabet[c]]
			p_x_t2 *= theta0[self.alphabet[c]]

		return (p_x_t1,p_x_t2)

	def _indicator(self,k,letter):	
		return 1 if self.alphabet[letter] == k else 0


	def _squash(self,z_array):
		"""
		Prob axioms. see Bailey and Elkan 1995
		"""
		array = z_array
		for i in range(self.n_seq):
			for j in range(len(self.sequenceList[i])):
				pass

		return array

	def _smooth(self,z_array):
		"""
		Prob axioms. see Bailey and Elkan 1995
		"""
		array = z_array
		for i in range(self.n_seq):
			for j in range(len(self.sequenceList[i])):
				pass

		return array
		
		
		
		
	def optimize_predict(self, n_iter = 10, beta = 0.01):
		"""
		Expectation Maximization Algorithm to find position weight matrix using all sequences in the file

		References:  Bailey, Elkan 1994. Link: https://www.ncbi.nlm.nih.gov/pubmed/7584402

		#############################################################################################

		Hyper Parameters:

		Number of iterations = 1000 (default)
		Dirichlet Distribution parameter (see eq 13) = beta = 0.01

		Process:

		Nucleic acids (DNA) : A, C, G, T	
		n = number of overlapping sequences formed from all sequences in the sequences.fa file. ~(seqLen - ML) * N

		Latent Variables =  lambda1 (Motif Class), lambda2 (Background Class) ---> Priors
							theta0 - array of length 4 (background model), 
							theta[1:ML] - array of shape MLx4  (motif model)
		
		Starting points: 
					Lambda1 - from Sqrt(n_seq)/n to 1/(2*ML)
					theta0 - computed from overall dataset - array of length 4
					theta1 - computed from the derived dataset (set of n overlapping sequences) - array of shape MLx4

		#############################################################################################

		Output:

		:predictedmotif
		:predictedsites - (one number for each of numSequnces lines)

		"""

		#############################################################################################
		## Initialization - search for starting points

		
		list_lambda = []
		avg_freq = [0]*self.lenAlphabet # Avg frequency of each of A,C,G,T in the dataset
		beta_arr = [] # beta * avg_freq

		for seq in self.sequenceList:
			for c in seq:
				avg_freq[self.alphabet[c]] += 1

		avg_freq = [i/np.sum(avg_freq) for i in avg_freq]
		beta_arr = [i * beta for i in avg_freq]

		# 2 Mixture Models
		overlappingSeq = [] # length = (seqLen - ML) * N
		for seq in self.sequenceList: # Motif class 
			for i in range(len(seq) - (self.motifLen - 1)):
				overlappingSeq.append(seq[i:i+self.motifLen]) 
		n = len(overlappingSeq)	

		
		from collections import Counter
		subseqctr = Counter()
		motif_site = []
		for o in range(len(overlappingSeq)):
			subseqctr[overlappingSeq[o]] += 1
			if overlappingSeq[o] == 'AGGATTA':
				motif_site.append(o % 494)
		top_m = subseqctr.most_common(5)
		dist_dict = Counter()
		for i in top_m:
			for j in top_m:
				if j[0] != i[0]:
					k = (i[0],j[0])
					for a,b in zip(i[0], j[0]):
						if a == b:
							dist_dict[k] += 1
		print(dist_dict)
		init_p = np.zeros(shape=(4, self.motifLen)) + 0.17
		for l in overlappingSeq:
			init_p = np.zeros(shape=(self.motifLen, 4)) + 0.17

			for j in range(self.motifLen):
				character = self.alphabet[l[j]]
				init_p[j][character] = 0.5


		theta0 = avg_freq.copy() # Background class; ex: [0.25,0.25,0.25,0.25] ; theta2 in paper1994
		theta1 = [[0]*self.lenAlphabet for l in range(self.motifLen)]		#theta1 in paper1994
		
		for seq in overlappingSeq:
			for i,c in enumerate(seq):
				theta1[i][self.alphabet[c]] +=1
		theta1 = [[ i/np.sum(pos) for i in pos] for pos in theta1] # Motif class ;similar to theta1.

		# print(theta0)
		# print(theta1)
		# Priors = [Lambda1 (motif class), lambda2 (background class)]; sum = 1
		
		prior_motif = 2e-3#;pass
		prior_background = 1-prior_motif#;pass

		# Assignments
		z_array = [[0]*(len(self.sequenceList[i]) - (self.motifLen - 1)) \
						for i in range(self.n_seq)] # Bayesian posterior; value = probability that sequence belongs to motif class. Ex: shape = 10x500
		
		# Erasing Factors
		e_array = [[1]*(len(self.sequenceList[i])) for i in range(self.n_seq)]

		"""Note: Z, z are same and E, e are same just expresseds in 2D,1D vectors respectively"""

		#############################################################################################
		## EM iterations
		
		while(n_iter>0):
			
			## Expectation step - (eq 4)
			for i in range(self.n_seq):
				for j in range(len(self.sequenceList[i])-self.motifLen):
					a, b = self._conditional_proba(overlappingSeq[10*i+j],theta0,theta1)
					z_array[i][j] = a*prior_motif/(a*prior_motif+b*prior_background) #eq4

			# normalizing z_array such that over each sequence i = 1 .. n  sum(z_i) <= 1. Probability axioms
			z_array = self._squash(z_array)
			z_array = self._smooth(z_array)

			# Update erasing factors
			for i in range(self.n_seq):
				for j in range(len(self.sequenceList[i])-self.motifLen):
					e_array[i][j] *= reduce(operator.mul, [1-z_array[i][k] for k in range(j-self.motifLen+1, j+1)], 1)
					
			## Maximization step - (eq 5 & 13)
			# priors
			prior_motif = np.sum(z_array)/n
			prior_background = 1- prior_motif

			# motifs
			c_back_arr = [0] * self.lenAlphabet # 1x4
			c_motif_arr = [[0] * self.lenAlphabet for l in range(self.motifLen)] # MLx4
			
			for k in range(self.lenAlphabet):
				for i in range(self.n_seq):
					for j in range(len(self.sequenceList[i])-self.motifLen):
						for l in range(self.motifLen):
							c_back_arr[k] += (1-z_array[i][j]) * self._indicator(k,overlappingSeq[10*i+j][l]) # eq9

				for l in range(self.motifLen):
					for i in range(self.n_seq):
						for j in range(len(self.sequenceList[i])-self.motifLen):
							c_motif_arr[l][k] += e_array[i][j]*z_array[i][j] * self._indicator(k,overlappingSeq[10*i+j][l])# eq10

			
			for j in range(self.lenAlphabet): # eq13
				theta0[j] = (c_back_arr[j] + beta_arr[j]) / (np.sum(c_back_arr) + beta)
				for i in range(len(theta1)):
					theta1[i][j] = (c_motif_arr[i][j] + beta_arr[j]) / (np.sum(c_motif_arr[i]) + beta) 

			n_iter-=1

		# print(theta0)
		# print(theta1)

		#############################################################################################
		## Found Motif and predicted Sites
		pred_sites = [0] * self.n_seq
		for i in range(self.n_seq):
			pred_sites[i] = z_array[i].index(max(z_array[i]))

		pred_motif = theta1
		#############################################################################################
		## return values
		return pred_sites, pred_motif


# """
# Tester code
# """	

dest_folder = 'results/dataset' + str(0)

sol = MotifFinder()
sol.set_motif_length(dest_folder + '/motiflength.txt')
sol.import_sequences(dest_folder + '/sequences.fa')
sites, motif = sol.optimize_predict()

print(sol.motifLen)
print(sol.sequenceList)
print(sol.n_seq)
print(sol.seqLen)
print(sites)
print(motif)

		 
