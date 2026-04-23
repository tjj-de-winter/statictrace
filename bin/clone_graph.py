# python code with graph networks

from sklearn.metrics import jaccard_score
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import itertools as it
import os,sys
from pandarallel import pandarallel
import numpy as np

G = nx.Graph()

csv = sys.argv[1]
jaccard = float(sys.argv[2]) # 0.3, 0.5

matrix = pd.read_csv(csv, index_col=0)

cells = matrix.index.to_list()

BC_dict = {cell: matrix.loc[cell, :].to_list() for cell in cells}

G.add_nodes_from(cells)

barcode_col = matrix.columns.to_list()

### Parameters ###
hamming_colapse = True
max_distance = 2

def distance(bc1, bc2):
	return np.sum([1 for b1, b2 in zip(bc1, bc2) if b1 != b2])

def hamming_check(bc_cell1, bc_cell2):
	barcodes1 = [bc for b, bc in zip(bc_cell1, barcode_col) if b == 1]
	barcodes2 = [bc for b, bc in zip(bc_cell2, barcode_col) if b == 1]

	collapse = []
	for i1, bc1 in enumerate(barcodes1):
	    for i2, bc2 in enumerate(barcodes2):
	        d = distance(bc1, bc2)
	        if d <= max_distance:
	            collapse.append((i1, i2))

	if len(collapse) > 0:
		for comb in collapse:
			bc1i, bc2i = barcode_col.index(barcodes1[comb[0]]), barcode_col.index(barcodes2[comb[1]])
			bc_cell1[bc1i] = 1
			bc_cell2[bc2i] = 1
	return bc_cell1, bc_cell2

def connection(comb, hamming_colapse=True):
	cell1, cell2 = comb
	bc_cell1  = BC_dict[cell1]
	bc_cell2  = BC_dict[cell2]
	
	if hamming_colapse:
		bc_cell1, bc_cell2 = hamming_check(bc_cell1, bc_cell2)

	JI = jaccard_score(bc_cell1, bc_cell2)
	if JI >= jaccard:
		return True
	else:
		return False

cell_combinations = list(it.combinations(cells, 2))
df = pd.DataFrame([cell_combinations],
				  index=['comb'],
				  columns=range(len(cell_combinations)),
				  ).T

pandarallel.initialize(progress_bar=True)
df['connection'] = df['comb'].parallel_apply(lambda comb: connection(comb))

df_edges = df[df['connection'] == True]

edges = df_edges['comb'].to_list()

G.add_edges_from(edges)

if hamming_colapse:
	csv = csv.replace('.csv', f'.hamming{max_distance}.csv')

outname = csv.replace('.csv',f'.graph.JI{jaccard}.gml')
nx.write_gml(G, outname)

print(f'###### FINAL EDGES = {G.number_of_edges()} ######')

fig, ax = plt.subplots(dpi=600)

isolates = list(nx.isolates(G))
non_isolates = [n for n in G.nodes if n not in isolates]

nx.draw(G, with_labels=False, node_size=8, ax=ax, nodelist=non_isolates)
outfilename = outname.replace('.gml','.png')
ax.set_title(outname.strip('.gml'))
fig.savefig(outfilename, bbox_inches='tight', dpi=600)

print(outfilename, 'made')
