### Description ###

# Given two fastqs containing lineage barcodes extracted from 10X scRNA-seq
# Get the cell barcode from R1 and the lineage barcode from R2

### Import packages ###

import sys, os
import argparse as argp
import re
import pandas as pd
import gzip
import numpy as np
import time
from itertools import combinations, product


### Input variables ###
parser = argp.ArgumentParser(description = 'Extract lineage barcode and cell barcode from paired-end fastq files')
parser.add_argument('--fq1', help = 'R1 fastq file, containing cell BC')
parser.add_argument('--fq2', help = 'R2 fastq file containing lineage BC')
parser.add_argument('--outdir', help = 'output directory', type = str, default = './')
parser.add_argument('--outprefix', help = 'outprefix')
args = parser.parse_args()

fq1 = args.fq1
fq2 = args.fq2
outprefix = args.outprefix
outdir = args.outdir

### Parameters ###
BC_1 = "GNNCTGNNCTNNACNNNNCGNNATNNGACNN" # barcode structure for BC1
BC_2 = "GNNGANNNGACNNTCNNGCNNACTNNNAGNN" # barcode structure for BC2

# BC_1 = "CAGNNCTGNNCTNNACNNNNCGNNATNNGACNN" # barcode structure for BC1
# BC_2 = "CAGNNGANNNGACNNTCNNGCNNACTNNNAGNN" # barcode structure for BC2

start = time.time()

### Functions ###

def reverse_complement(sequence):
    complement = {'A':'T',
                  'T':'A',
                  'C':'G',
                  'G':'C',
                  'N':'N'}
    return ''.join([complement[n] for n in sequence[::-1]])

# def hamming_distance(BC_seq):
#     '''Generate all possible barcode sequences with 1 hamming distance difference'''
#     nucleotides = ['A', 'C', 'G', 'T']
#     neighbors = [BC_seq]

#     for i, n in enumerate(BC_seq):
#         for alt in nucleotides:
#             if alt != n and alt != 'N':
#                 neighbor = BC_seq[:i] + alt + BC_seq[i+1:]
#                 neighbors.append(neighbor)
#     return neighbors

from itertools import combinations, product

def hamming_distance(seq, max_dist=2):
    """
    Generate all sequences within Hamming distance <= max_dist
    of the given sequence. 'N' positions are skipped (treated as wildcards).
    """
    nucleotides = ['A', 'C', 'G', 'T']
    neighbors = {seq}  # include original

    L = len(seq)

    # Distance = 1
    for i in range(L):
        if seq[i] == 'N':
            continue
        for nc in nucleotides:
            if nc != seq[i]:
                neighbors.add(seq[:i] + nc + seq[i+1:])

    # Distance = 2
    if max_dist >= 2:
        for i, j in combinations(range(L), 2):
            if seq[i] == 'N' or seq[j] == 'N':
                continue
            for nc1, nc2 in product(nucleotides, repeat=2):
                if nc1 == seq[i] or nc2 == seq[j]:
                    continue
                neighbors.add(seq[:i] + nc1 + seq[i+1:j] + nc2 + seq[j+1:])

    return neighbors

def detect_barcode(seq):

    indices = indices_BC1

    for pattern in patterns_BC1:
        match = pattern.search(seq)
        if match:
            barcode = ''.join([match[0][i] for i in indices])
            index_start = match.start()
            return 'BC1', barcode, index_start

    indices = indices_BC2

    for pattern in patterns_BC2:
        match = pattern.search(seq)
        if match:
            barcode = ''.join([match[0][i] for i in indices])
            index_start = match.start()
            return 'BC2', barcode, index_start
    return 'no BC detected', '', 0

def get_phred_scores(index_start, phred_score_read, bc_type):
    indices_dict = {'BC1': indices_BC1, 'BC2': indices_BC2}
    indices = indices_dict[bc_type]
    phred_score_read_rev = phred_score_read[::-1]
    phred_score_read_rev = phred_score_read_rev[index_start:]
    barcode_phred = ''.join([phred_score_read_rev[i] for i in indices])
    return barcode_phred

### Code ###

indices_BC1 = [i for i, x in enumerate(BC_1) if x == "N"]
indices_BC2 = [i for i, x in enumerate(BC_2) if x == "N"]

BC1s = hamming_distance(BC_1)
patterns_BC1 = [re.compile(bc.replace('N','[ACGT]')) for bc in BC1s]

BC2s = hamming_distance(BC_2)
patterns_BC2 = [re.compile(bc.replace('N','[ACGT]')) for bc in BC2s]

print('len BC1s:', len(BC1s), 'len BC2s:', len(BC2s))

# Create output directory if it does not exist
if not os.path.isdir(outdir):
    os.system('mkdir '+outdir)

n_bc1 = 0
n_bc2 = 0

f1 = gzip.open(fq1)
f2 = gzip.open(fq2)

outfile = os.path.join(outdir, f'{outprefix}.linBCs.csv')
n_reads = 0

if os.path.exists(outfile):
    os.remove(outfile)

with open(outfile, 'a', buffering=1) as fout:
    fout.write(f'readname,cellBC,cellBC_phred,linBC,linBC_phred,BC_structure\n')
    for idx, (l1, l2) in enumerate(zip(f1.readlines(), f2.readlines())):
        try:
            l1 = str(l1.rstrip().rsplit()[0], 'utf-8')
            l2 = str(l2.rstrip().rsplit()[0], 'utf-8')
        except:
            continue
        l = np.mod(idx,4)
        if l1.startswith('@'):
            readname = l2.replace(' ','')
        elif l == 1:
            n_reads += 1
            bc_type, bc, index_start = detect_barcode(reverse_complement(l2))
            cellbc = l1
        elif l == 3:
            if bc_type != 'no BC detected':
                bc_phred = get_phred_scores(index_start, l2, bc_type)
                fout.write(f'{readname},{cellbc},{l1},{bc},{bc_phred},{bc_type}\n')
                # print(f'barcode detected: {l1} - {l2}')
                if bc_type == 'BC1':
                    n_bc1 += 1
                elif bc_type == 'BC2':
                    n_bc2 += 1
        fout.flush()

print(time.time()-start)

print(f'processing {fq1} + {fq2} done')
print(f'Total number of reads processed: {n_reads}')
print(f'Total number of barcodes detected: {n_bc1 + n_bc2}')
print(f'Percentage of barcodes detected: {round((n_bc1 + n_bc2)/n_reads*100, 2)}%')
print(f'Number of barcodes with BC1 structure: {n_bc1}')
print(f'Number of barcodes with BC2 structure: {n_bc2}')
print(f'{outfile} generated')
