### Description ###

# add clonal family information to an adata object and save as h5ad file

### Import packages ###

import scanpy as sc
import pandas as pd
import numpy as np
import argparse as argp


### Input variables ###
parser = argp.ArgumentParser(description = 'add clonal family information to an adata object and save as h5ad file')
parser.add_argument('--h5ad', help = 'h5ad file without clonal information') #10X/h5ad/PSC694_v1.h5ad
parser.add_argument('--h5ad_out', help = 'output h5ad file with clonal information') #10X/h5ad/PSC694_clonal_v1.h5ad'
parser.add_argument('--clones_BC1', help = '.BC1.clones.csv file') # v7.linBC_matrix.BC1.clones.csv'
parser.add_argument('--clones_BC2', help = '.BC2.clones.csv file') # v7.linBC_matrix.BC2.clones.csv

args = parser.parse_args()

h5ad_in = args.h5ad
h5ad_out = args.h5ad_out
clones_BC1 = args.clones_BC1
clones_BC2 = args.clones_BC2

### Code ###

adata = sc.read_h5ad(h5ad_in)

df_bc1 = pd.read_csv(clones_BC1, names=['cellbc', 'clonal_family'], index_col=0)
df_bc2 = pd.read_csv(clones_BC2, names=['cellbc', 'clonal_family'], index_col=0)

def get_clone(cellbc, df_bc):
	bc_dict = {cellbc: clone for cellbc, clone in zip(df_bc.index, df_bc.clonal_family)}

	cellbc = cellbc.strip('-1')
	if cellbc in bc_dict.values():
		clone = bc_dict[cellbc]
		print(clone)
	else:
		clone = np.nan
	return clone

cellbcs = adata.obs.index
adata.obs['BC1_clone'] = [get_clone(cellbc, df_bc1) for cellbc in cellbcs]
adata.obs['BC2_clone'] = [get_clone(cellbc, df_bc2) for cellbc in cellbcs]

adata.obs['BC1_clone'] = adata.obs['BC1_clone']
adata.obs['BC2_clone'] = adata.obs['BC2_clone']

adata.write(h5ad_out)

print(len(set(df_bc1['clonal_family'])),'clones added for BC1')
print(len(set(df_bc2['clonal_family'])),'clones added for BC2s')
print('data can be found in:', h5ad_out)
