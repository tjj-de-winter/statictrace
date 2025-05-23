### Description ###

# Add linBC consensus sequences to a matrix (cellBCs vs linBCs) and export as CSV

### Import packages ###

import sys, os
import pandas as pd
import glob
import numpy as np

### Input variables ###
parser = argp.ArgumentParser(description = 'Add linBC consensus sequences to a matrix (cellBCs vs linBCs) and export as CSV')
parser.add_argument('--starcode_dir', help = 'folder containing .true_barcodes.txt files')
parser.add_argument('--outdir', help = 'output directory', type = str, default = './')
parser.add_argument('--outprefix', help = 'outprefix')
args = parser.parse_args()

path = args.starcode_dir
fq2 = args.fq2
outprefix = args.outprefix
outdir = args.outdir

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

### Code ###

BC1_files = glob.glob(path+'/*.BC1.true_barcodes.txt')
BC2_files = glob.glob(path+'/*.BC2.true_barcodes.txt')

df_bc1 = generate_matrix(BC1_files)

df_bc1.to_csv(f'{outdir}{outprefix}.BC1.linBC_matrix.csv')

df_bc2 = generate_matrix(BC2_files)

df_bc2.to_csv(f'{outdir}{outprefix}.BC2.linBC_matrix.csv')