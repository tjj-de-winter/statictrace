# statictrace

This pipeline assumes cells are barcoded using two rounds of lentivirus based barcoding (i.e. BC1 for round 1 and BC2 for round 2), cells are processed using 10X scRNA-seq.

<b>Step 1</b>: Perform QC and process the scRNAseq dataset in scanpy, and export the adata object as an h5ad file

<b>Step 2</b>: Extract all lineage barcodes and corresponding cell barcodes from a FASTQ file (first run trimgalore to remove low quality reads and bases), a mismatch of 1 nucleotide in the constant nucleotide part of the lineage barcode is allowed. After extracting the barcode, the constant nucleotides are removed and only random nucleotides are kept, resulting in a 16 basepair lineage barcode. This script generates an output CSV file called <i>.linBCs.csv</i>

```
python extract_BCs.v2fast.py --fq1 <fastq R1 file> --fq2 <fastq R1 file> --outdir <output path> --outprefix <output prefix>
```

<b>Step 3</b>: Get all consensus sequences for lineage tracing barcodes extracted in step 2 that are present in the cells that passed QC (step 1), a mismatch of 1 nuleotide in the cell barcode is allowed. Default settings for allowable levensthein distance to compute consensus sequence is set to 3. This script generates four output files (two for BC1 and two for BC2) for each cell barcode.

```
python consensus_linBC.py --CSVdir <path to .linBCs.csv> --h5ad <h5ad file> --distance <levensthein distance> --outdir <output path>
```

<b>Step 4</b>: Extract true barcodes based on starcode consensus read clusters. Default minimum read depth to call a lineage barcode is set to 4 reads. This script generates a <i>.BC1.true_barcodes.txt</i> file and <i>.BC2.true_barcodes.txt</i> file for each cell, note that the file is empty if no true lineage barcode could be detected.
```
./consensus_linBC.sh <path to .clustered.txt> <n_reads>
```

<b>Step 5</b>: Generate a matrix containing cells (cell barcodes) as index and true lineage barcodes as columns. This script generates <i>.BC1.linBC_matrix.csv</i> and <i>.BC2.linBC_matrix.csv</i> files.
```
python linBC_matrix.py --starcode_dir <path to .true_barcodes.txt> --outdir <output path> --outprefix <output prefix>
```

<b>Step 6</b>: Group cells in clonal families by barcode presence using jaccard similarity scoring, cluster the cells using hierarchical clustering. The detault minimum jaccard filtering parameter for clustering is set to 0.5. The minimum clone size is set to 2 cells per clone. This script generates a <i>.clones.csv</i> file
```
python clonal_families --BCmatrix <.linBC_matrix.csv> --outdir <output path>
```

<b>Step 7</b>: Append clonal information from step 6 to the adata h5ad file generated in step 1
```
python appened_clones_h5ad.py --h5ad <h5ad file> --h5ad_out <h5ad output file> --clones_BC1 <.BC1.clones.csv> --clones_BC2 <.BC2.clones.csv>
```



