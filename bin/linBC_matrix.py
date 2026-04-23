### Description ###

# Add linBC consensus sequences to a matrix (cellBCs vs linBCs) and export as CSV

### Import packages ###

import sys, os
import pandas as pd
import glob
import numpy as np
import argparse as argp

### Input variables ###
parser = argp.ArgumentParser(description = 'Add linBC consensus sequences to a matrix (cellBCs vs linBCs) and export as CSV')
parser.add_argument('--starcode_dir', help = 'folder containing .true_barcodes.txt files')
parser.add_argument('--outdir', help = 'output directory', type = str, default = './')
parser.add_argument('--outprefix', help = 'outprefix')
args = parser.parse_args()

path = args.starcode_dir
outprefix = args.outprefix
outdir = args.outdir

if not os.path.exists(outdir):
	os.mkdir(outdir)

### Parameters ###
min_BCs = 1 # minimum length of barcode set
max_BCs = 10 # maximum length of barcode set

### Functions ###

def generate_matrix(files):

	barcode_dict = {}
	all_linBCs = []
	for file in files:
		txt = open(file, 'r')
		linBCs = []
		cellBC = os.path.basename(file).split('.')[0]
		for line in txt.readlines():
			line = line.strip('\n')
			linBCs.append(line)
			all_linBCs.append(line)
		if len(linBCs) > 0:
			barcode_dict[cellBC] = linBCs

	index = list(barcode_dict.keys())
	columns = list(set(all_linBCs))
	zeros = np.zeros((len(index),len(columns)))
	df = pd.DataFrame(zeros, index=index, columns=columns)

	for cellBC in barcode_dict.keys():
		linBCs = barcode_dict[cellBC]
		for linBC in linBCs:
			df.loc[cellBC, linBC] = 1

	return df

def filter(df):
	return df[(df.sum(axis=1)>=min_BCs)&(df.sum(axis=1)<=max_BCs)]

### Code ###

BC1_files = glob.glob(path+'/*.BC1.true_barcodes.txt')
BC2_files = glob.glob(path+'/*.BC2.true_barcodes.txt')



df_bc1 = generate_matrix(BC1_files)
df_bc1 = filter(df_bc1)
df_bc1.to_csv(os.path.join(outdir,f'{outprefix}.BC1.linBC_matrix.csv'))

df_bc2 = generate_matrix(BC2_files)
df_bc2 = filter(df_bc2)
df_bc2.to_csv(os.path.join(outdir,f'{outprefix}.BC2.linBC_matrix.csv'))
