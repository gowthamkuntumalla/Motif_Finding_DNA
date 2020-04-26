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

class MotifFinder:
	""" Solution Class """
	
	def __init__(self):
		self.sequenceList = []
		self.motifLen = 0
	
	
	def set_motif_length(self,lengthfile):
		""" 
		number of DNA sequences in the file sequences.fa 
		"""
		with open(lengthfile,'r') as f:
			lines = f.readlines()
		self.motifLen = lines[0]
	
	
	def import_sequences(self,sequencesFile):
		""" 
		sequencesFile is the file name in results/dataset{i}
		"""
		with open(sequencesFile,'r') as f:
			lines = f.readlines()

		for i,line in enumerate(lines):
			if i%2 == 1:
				self.sequenceList.append(line.strip())

	
	def predict_motif(self):
		"""
		Algorithm to find position weight matrix using all sequences in the file
		output
		:predictedmotif.txt
		"""
		pass

	
	def predict_sites(self):
		"""
		uses predictedmotif to locate one site per sequence

		output
		:predictedsites.txt - (one number for each of numSequnces lines)	
		"""
		pass


"""
Tester code
"""	

# dest_folder = 'results/dataset' + str(0)

# sol = MotifFinder()
# sol.set_motif_length(dest_folder + '/motiflength.txt')
# sol.import_sequences(dest_folder + '/sequences.fa')
# sol.predict_motif() # creates predictedmotif.txt
# sol.predict_sites() # creates predictedsites.txt

# print(sol.motifLen)
# print(sol.sequenceList)
		
		 
