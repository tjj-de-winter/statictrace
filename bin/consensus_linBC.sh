#!/bin/bash

# Get true barcodes from starcode generated cluster files, only select consensus barcodes with a minimal read count 

folder=$1 #starcode
n_reads=4 # minimum read depth to call a true barcode

cd $folder

for file in $(ls *.clustered.txt)
do
	outfile=${file%.clustered.txt}.true_barcodes.txt
	awk -v threshold="$n_reads" '$2 > threshold {print $1}' "$file" > $outfile
done