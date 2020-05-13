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
import copy
import os
import numpy as np
from functools import reduce
import operator
from collections import Counter


class MotifFinder:
    """ Solution Class """

    def __init__(self):
        self.sequenceList = []
        self.motifLen = 0
        self.n_seq = 0  # number of sequences in the file sequences.fa
        self.seqLen = 0  # avg length of a sequence.
        self.alphabet = {'A': 0, 'G': 1, 'C': 2, 'T': 3}  # positions of Nucelic Acid alphabet
        self.lenAlphabet = len(self.alphabet.keys())

    def set_motif_length(self, lengthfile):
        """
        number of DNA sequences in the file sequences.fa
        """
        with open(lengthfile, 'r') as f:
            lines = f.readlines()
        self.motifLen = int(lines[0])

    def import_sequences(self, sequencesFile):
        """
        sequencesFile is the file name in results/dataset{i}
        """
        with open(sequencesFile, 'r') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if i % 2 == 1:
                self.sequenceList.append(line.strip())

        self.n_seq = len(self.sequenceList)
        self.seqLen = int(np.mean([len(l) for l in self.sequenceList]))

    def _conditional_proba(self, seq, theta0, theta1):
        """
        Private func.
        To be called in expectation step of iterations

        output:
            tuple(P(Xi|theta1),P(Xi|theta2))
        """

        p_x_t1 = 1  # motif class
        p_x_t2 = 1  # background class

        for i, c in enumerate(seq):
            p_x_t1 *= theta1[i][self.alphabet[c]]
            p_x_t2 *= theta0[self.alphabet[c]]

        return (p_x_t1, p_x_t2)

    def _indicator(self, k, letter):
        return 1 if self.alphabet[letter] == k else 0

    def _squash(self, z_array):
        """
        Prob axioms. see Bailey and Elkan 1995
        """
        array = z_array
        for i in range(self.n_seq):
            for j in range(len(self.sequenceList[i])):
                pass

        return array

    def _smooth(self, z_array):
        """
        Prob axioms. see Bailey and Elkan 1995
        """
        array = z_array
        for i in range(self.n_seq):
            for j in range(len(self.sequenceList[i])):
                pass

        return array

    ## See Appendix From Bailey-Elkan 1995
    ## Calculates log likelihood than normalizes to ONE
    def calculate_z_table(self, p, overlapseqs):

        z_matrix = np.zeros((self.n_seq, len(self.sequenceList[0]) - self.motifLen + 1))
        for i in range(len(overlapseqs)):
            prob = 1
            for j in range(self.motifLen):
                prob *= p[j][self.alphabet[overlapseqs[i][j]]]
            z_matrix[int(i / len(z_matrix[0]))][i % len(z_matrix[0])] = prob
        z_matrix = z_matrix / z_matrix.sum(axis=1, keepdims=1)
        return z_matrix

    def optimize_predict(self, n_iter=100, beta=0.05):
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
        avg_freq = [0] * self.lenAlphabet  # Avg frequency of each of A,C,G,T in the dataset
        beta_arr = []  # beta * avg_freq

        for seq in self.sequenceList:
            for c in seq:
                avg_freq[self.alphabet[c]] += 1

        avg_freq = [i / np.sum(avg_freq) for i in avg_freq]
        beta_arr = [i * beta for i in avg_freq]

        # 2 Mixture Models
        overlappingSeq = []  # length = (seqLen - ML) * N
        seqc = 0
        for seq in self.sequenceList:  # Motif class
            for i in range(len(seq) - (self.motifLen - 1)):
                overlappingSeq.append(seq[i:i + self.motifLen])
            seqc += 1
        n = len(overlappingSeq)

        """
        Find a good P starting Point, Loop through the most common Subsequences and maximize log likelihood
        The motif will Probably be in multiple Sequences. 
        Description of how to initialize is in Bailey Elkan 1995, section 2.1
        """

        subseqctr = Counter()
        for o in range(len(overlappingSeq)):
            subseqctr[overlappingSeq[o]] += 1
        top_m = subseqctr.most_common(10) # top 10 motif candidates

        max_likely = 0
        best_pwm = np.zeros(shape=(self.motifLen, self.lenAlphabet)) + 0.17
        best_motif = None
        best_z_array = None
        for candidate_motif in top_m:
            candidate_start = np.zeros(shape=(self.motifLen, self.lenAlphabet)) + 0.17

            for j in range(self.motifLen):
                character = self.alphabet[candidate_motif[0][j]]
                candidate_start[j][character] = 0.5
            cand_z_array = self.calculate_z_table(candidate_start, overlappingSeq)
            prior_motif = np.sum(cand_z_array) / n
            max_likelihood = np.sum(np.max(cand_z_array, axis=1))
            if max_likelihood > max_likely:
                max_likely = max_likelihood
                best_pwm = candidate_start
                best_motif = candidate_motif
                best_z_array = cand_z_array

        """
        Alternatively, loop through all subsequences. 
        Commented out, not in use. 
        """
        # for candidate_motif in overlappingSeq:
        #     candidate_start = np.zeros(shape=(self.motifLen, 4)) + 0.17
        #
        #     for j in range(self.motifLen):
        #         character = self.alphabet[candidate_motif[j]]
        #         candidate_start[j][character] = 0.5
        #     cand_z_array = self.calculate_z_table(candidate_start, overlappingSeq)
        #     prior_motif = np.sum(cand_z_array) / n
        #     max_likelihood = np.sum(np.max(cand_z_array, axis=1))
        #     if max_likelihood > max_likely:
        #         max_likely = max_likelihood
        #         best_pwm = candidate_start
        #         best_motif = candidate_motif
        #     print('hi')

        # theta0 = avg_freq.copy()  # Background class; ex: [0.25,0.25,0.25,0.25] ; theta2 in paper1994
        '''
        Copy over the "estimated" best Motif Matrix from above to theta. 
        
        Really should do a MAXIMIZATION Step first, and then check max likelihood? 
        
        '''
        theta1 = copy.deepcopy(best_pwm)  # theta1 in paper1994, motif_class

        # print(theta0)
        # print(theta1)
        # Priors = [Lambda1 (motif class), lambda2 (background class)]; sum = 1

        prior_motif = 2e-3  # ;pass
        prior_background = 1 - prior_motif  # ;pass

        # Assignments
        z_array = best_z_array  # Bayesian posterior; value = probability that sequence belongs to motif class. Ex: shape = 10x500

        # Erasing Factors
        # e_array = [[1] * (len(self.sequenceList[i])) for i in range(self.n_seq)]

        """Note: Z, z are same and E, e are same just expresseds in 2D,1D vectors respectively"""

        #############################################################################################
        ## EM iterations

        while (n_iter > 0):

            ## Expectation step - (eq 4)
            '''
            Using the simplified E step formula from 1995 Appendix. Thus, lambdas not being used
            Is this a possible improvement?
            '''

            # normalizing z_array such that over each sequence i = 1 .. n  sum(z_i) <= 1. Probability axioms
            ### Not using  ####
            # z_array = self._squash(z_array)
            # z_array = self._smooth(z_array)

            '''
            Erasing Was mentioned in the 1994 paper, but it seems to be for finding multiple Motifs ?
            i.e start erasing a motif so that new ones can arise 
            Not using right now,
            Update erasing factors - Do we need this? Seems to be for multiple Motifs.
            '''
            # for i in range(self.n_seq):
            #     for j in range(len(self.sequenceList[i]) - self.motifLen):
            #         e_array[i][j] *= reduce(operator.mul,
            #                                 [1 - z_array[i][k] for k in range(j - self.motifLen + 1, j + 1)], 1)

            ## Maximization step - (eq 5 & 13)
            # priors
            prior_motif = np.sum(z_array) / n
            prior_background = 1 - prior_motif

            # motifs
            c_back_arr = [0] * self.lenAlphabet  # 1x4
            c_motif_arr = [[0] * self.lenAlphabet for l in range(self.motifLen)]  # MLx4

            for k in range(self.lenAlphabet):

                '''
                Not sure if we need to update the background formula, we know (i think?) that this is drawn
                from a random distribution
                '''
                # for i in range(self.n_seq):
                #     for j in range(len(self.sequenceList[i]) - self.motifLen):
                #         for l in range(self.motifLen):
                #             c_back_arr[k] += (1 - z_array[i][j]) * self._indicator(k,
                #                                                                    overlappingSeq[10 * i + j][l])  # eq9

                for l in range(self.motifLen):
                    for i in range(self.n_seq):
                        for j in range(len(self.sequenceList[i]) - self.motifLen + 1):
                            c_motif_arr[l][k] += z_array[i][j] * \
                                                 self._indicator(k, overlappingSeq[((len(self.sequenceList[i]) - self.motifLen + 1) * i) + j][l])  # eq10
                # for l in range(self.motifLen):
                #     for i, j in zip(range(self.n_seq), curr_pos_guess):
                #         c_motif_arr[l][k] += z_array[i][j] * self._indicator(k, overlappingSeq[(i * 494) + j][l])

            # for j in range(self.lenAlphabet):  # eq13
            #             #     theta0[j] = (c_back_arr[j] + beta_arr[j]) / (np.sum(c_back_arr) + beta)
            #             #     for i in range(len(theta1)):
            #             #         theta1[i][j] = (c_motif_arr[i][j] + beta_arr[j]) / (np.sum(c_motif_arr[i]) + beta)
            '''
            Update the PWM
            '''
            for i in range(len(theta1)):
                for j in range(self.lenAlphabet):
                    theta1[i][j] = (c_motif_arr[i][j] + beta_arr[j]) / (np.sum(c_motif_arr[i]) + beta) # just so that 0 prob doesn't cause computational issues

            z_array = self.calculate_z_table(theta1, overlappingSeq)   #eq4

            n_iter -= 1

        # print(theta0)
        # print(theta1)

        #############################################################################################
        ## Found Motif and predicted Sites
        pred_sites = [0] * self.n_seq
        test = self.calculate_z_table(best_pwm, overlappingSeq)

        for i in range(self.n_seq):
            pred_sites[i] = np.where(z_array[i] == np.max(z_array[i]))[0][0]
        print(pred_sites)
        # for i in range(self.n_seq):
        #     pred_sites[i] = np.where(test[i] == np.max(test[i]))[0][0] # Is there change from the original guess?
        pred_motif = theta1
        print(pred_motif)
        #############################################################################################
        ## return values
        return pred_sites, pred_motif


# """
# Tester code
# """

# dest_folder = 'results/dataset' + str(12)

# sol = MotifFinder()
# sol.set_motif_length(dest_folder + '/motiflength.txt')
# sol.import_sequences(dest_folder + '/sequences.fa')
# sites, motif = sol.optimize_predict()

# print(sol.motifLen)
# print(sol.sequenceList)
# print(sol.n_seq)
# print(sol.seqLen)
# print(sites)
# print(motif)
