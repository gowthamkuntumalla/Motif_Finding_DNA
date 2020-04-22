""" 
Step1: Benchmarking 

Creating datasets for testing motif finding algorithm

Expected arguments:

	(a) A positive number called ICPC (“information content per column”) 
	(b) A positive integer called ML (“motif length”)
	(c) A positive integer called SL (“sequence length”)
	(d) A positive integer called SC (“sequence count”)

"""

import numpy as np
import scipy.optimize as opt
import os
import math

info_cont_per_col = [1.0, 1.5, 2.0]
motif_length = [6, 7, 8]
sequence_length = 500
sequence_count = [5, 10, 20]
nucleotides = ['A', 'G', 'C', 'T']
default_motif_len = 7
default_seq_ct = 10
default_icpc = 1.5


def gen_rand_seq(seq_len):
    seed = np.random.uniform(0, 4, seq_len)
    seq = []
    for i in range(seq_len):
        seq.append(nucleotides[int(seed[i])])
    return np.array(seq)


"""
Basically split the summation function in half, then solve for the root. 
"""
def infocontentcalc(x, dyn_sum, icpc, p1, p2):
    out = (x * np.log2(4 * x)) + ((1 - p1 - p2 - x) * np.log2(4 * (1 - p1 - p2 - x))) + dyn_sum - icpc
    return out


def generate_pwm(ml, icpc):
    motif = []
    if icpc == 2:
        for col in range(ml):
            winner = np.random.randint(0, 4)
            column = [0, 0, 0, 0]
            column[winner] = 1
            motif.append(column)
        return np.array(motif)

    for col in range(ml):
        p1, p2, p3, p4 = -1, -1, -1, -1
        err = .02
        calculated_icpc = icpc - (2 * err)
        while any(x < 0 for x in [p1, p2, p3, p4]) or calculated_icpc < icpc - err or calculated_icpc > icpc + err:
            p1 = np.random.uniform(.001, .06)
            p2 = np.random.uniform(.001, .06)
            p34 = 1 - p2 - p1
            sum_p12 = (p1 * np.log2(4 * p1)) + (p2 * np.log2(4 * p2))
            roots = opt.fsolve(infocontentcalc, .9, (sum_p12, icpc, p1, p2))
            p3 = 1 - roots[0] - p1 - p2
            p4 = roots[0]
            calculated_icpc = ic(p1, p2, p3, p4)
        column = np.array([p1, p2, p3, p4])
        np.random.shuffle(column)
        motif.append(column)
    return np.array(motif)


"""
Just a test function to make sure i had the information formula correct
"""
def ic(q, r, s, t):
    return (q * np.log2(4 * q)) + (r * np.log2(4 * r)) + (s * np.log2(4 * s)) + (t * np.log2(4 * t))


def generate_motif_str(ml, pwm):
    mstr = []
    for i in range(ml):
        mstr.append(np.random.choice(nucleotides, p=pwm[i]))
    return np.array(mstr)


def plant_in_string(ml, motifstr, sample):
    site = np.random.randint(0, len(sample) - ml)
    sample[site:site + ml] = motifstr
    return sample, site


os.mkdir('benchmarks')


def write_dataset_to_file(dataset_ctr, samples, pwm, sites, motif_len):
    os.mkdir('benchmarks/dataset' + str(dataset_ctr))
    # dataset_ctr += 1
    f = open('benchmarks/dataset' + str(dataset_ctr) + '/sequences.fa', 'w+')
    for idx, seq in enumerate(samples):
        f.write('>sequence' + str(idx) + '\n')
        out = ''
        for n in seq:
            out += n
        f.write(out + '\n')
    f.close()
    f = open('benchmarks/dataset' + str(dataset_ctr) + '/sites.txt', 'w+')
    for site_num in sites:
        f.write(str(site_num) + '\n')
    f.close()
    f = open('benchmarks/dataset' + str(dataset_ctr) + '/motif.txt', 'w+')
    f.write('>MOTIF{}    {}\n'.format(str(dataset_ctr), str(motif_len)))

    for col in pwm:
        col_str = ' '.join([str(entry) for entry in col])
        f.write(col_str + '\n')
    f.write('<')
    f.close()
    f = open('benchmarks/dataset' + str(dataset_ctr) + '/motiflength.txt', 'w+')
    f.write(str(motif_len))
    f.close()


dataset_ctr = 0

for x in info_cont_per_col:
    for i in range(10):
        samples = []
        motif_strs = []
        sites = []
        pwm = generate_pwm(default_motif_len, x)
        for s in range(default_seq_ct):
            samples.append(gen_rand_seq(sequence_length))
            motif_strs.append(generate_motif_str(default_motif_len, pwm))
        for strnum in range(len(samples)):
            samples[strnum], site = plant_in_string(default_motif_len, motif_strs[strnum], samples[strnum])
            sites.append(site)
        write_dataset_to_file(dataset_ctr, samples, pwm, sites, default_motif_len)
        dataset_ctr += 1

for m_len in [6, 8]:
    for i in range(10):
        samples = []
        motif_strs = []
        sites = []
        pwm = generate_pwm(m_len, default_icpc)
        for s in range(default_seq_ct):
            samples.append(gen_rand_seq(sequence_length))
            motif_strs.append(generate_motif_str(m_len, pwm))
        for strnum in range(len(samples)):
            samples[strnum], site = plant_in_string(m_len, motif_strs[strnum], samples[strnum])
            sites.append(site)
        write_dataset_to_file(dataset_ctr, samples, pwm, sites, m_len)
        dataset_ctr += 1

for seq_ct in [5, 20]:
    for i in range(10):
        samples = []
        motif_strs = []
        sites = []
        pwm = generate_pwm(default_motif_len, default_icpc)
        for s in range(seq_ct):
            samples.append(gen_rand_seq(sequence_length))
            motif_strs.append(generate_motif_str(default_motif_len, pwm))
        for strnum in range(len(samples)):
            samples[strnum], site = plant_in_string(default_motif_len, motif_strs[strnum], samples[strnum])
            sites.append(site)
        write_dataset_to_file(dataset_ctr, samples, pwm, sites, default_motif_len)
        dataset_ctr += 1
