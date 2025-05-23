### Description ###

# find clones by computing the jaccard index between sets of lineage barcodes

### Import packages ###

import pandas as pd
from sklearn.metrics import jaccard_score
import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
import argparse as argp

### Input variables ###
parser = argp.ArgumentParser(description = 'Extract clonal families using barcode jaccard similarity and hierarchical clustering')
parser.add_argument('--BCmatrix', help = 'BC matrix generated with python linBC_matrix.py')
parser.add_argument('--outdir', help = 'output directory', type = str, default = './')
args = parser.parse_args()

file = args.BCmatrix
outdir = args.outdir

### Parameters ###

jaccard_treshold = 0.5
min_clonesize = 2

### Functions ###

def jaccard_similarity(df):
	cells = list(df.index)
	df_jaccard = pd.DataFrame(index=cells,columns=cells).astype(float)
	combinations = []
	for cellbc1 in cells:
		for cellbc2 in cells:
			if tuple(sorted([cellbc1,cellbc2])) not in combinations:
				combinations.append(tuple(sorted([cellbc1,cellbc2])))
				barcode_list1 = df.loc[cellbc1]
				barcode_list2 = df.loc[cellbc2]
				jaccard = jaccard_score(barcode_list1, barcode_list2)
				df_jaccard.loc[cellbc1,cellbc2] = jaccard
				df_jaccard.loc[cellbc2,cellbc1] = jaccard
	return df_jaccard

def cluster(df, method='ward', optimal_ordering=True):
    data_array = df.to_numpy()
    
    dend_row = shc.dendrogram(shc.linkage(data_array, 
                                          method=method,
                                          optimal_ordering=optimal_ordering
                                         ), get_leaves=True, no_plot=True)

    df_clustered = df.iloc[dend_row['leaves'],dend_row['leaves']].copy()
    
    fig, ax = plt.subplots(dpi=300)

    im = ax.imshow(df_clustered, cmap='Greens')
    plt.colorbar(im, label='jaccard similarity')
    ax.set_title(f'barcode sharing')
    # plt.show()
    fig.savefig(file.replace('.csv','.jaccard.png'), bbox_inches='tight', dpi=300)
    plt.close(fig)
    return df_clustered

def group_clonal_families(df_clustered):
    processed_cells = set()
    clone_dict = {}
    n = 1
    for i, cell in enumerate(df_clustered.columns):
        if cell not in processed_cells:
            cloneID = f'{BC_structure}_clone{n}'
            clonal_cells = list(df_clustered[df_clustered[cell]>0].index)
            if len(clonal_cells) >= min_clonesize:
                clone_dict[cloneID] = clonal_cells
                n += 1
            processed_cells = processed_cells.union(set(clonal_cells))
    
    return clone_dict

def clonedict_to_csv(outfile):
    with open(outfile, 'w') as fout:
        for clone in clone_dict:
            cells = clone_dict[clone]
            for cell in cells:
                fout.write(f'{cell},{clone}\n')

### Code ###

BC_structure = file.split('.')[2]

df = pd.read_csv(file, index_col=0)
df = df.loc[:,df.sum(axis=0) > 1].copy() # remove orphan lineage barcodes
df = df.loc[df.sum(axis=1) > 0,:].copy() # remove cells that have no more barcodes

print('start calculating jaccard similarity')
df_jaccard = jaccard_similarity(df)
print('done calculating jaccard similarity')

df_jaccard.to_csv(file.replace('.csv','.jaccard.csv'))

# df_jaccard = pd.read_csv(file.replace('.csv','.jaccard.csv'), index_col=0)

df_jaccard[df_jaccard <= jaccard_treshold] = 0
print('start clustering data')
df_clustered = cluster(df_jaccard, method='ward', optimal_ordering=True)
clone_dict = group_clonal_families(df_clustered)
print('start grouping data')
clonedict_to_csv({outdir}+file.replace('.csv','.clones.csv'))
