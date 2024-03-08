#!/bin/bash

# This script organizes the RNA matrix, features and barcodes
# that were downloaded from NIH GEO GSE212447
# Into the correct directories so they can be read in with Seurat


cd data/rna_GSE212447/GSE212447_RAW/ || exit

samples=$(ls *.mtx.gz | sed 's/_matrix\.mtx\.gz$//'  | sort -u)

for sample in $samples; do
	mkdir -p "$sample"
	
	mv "${sample}_matrix.mtx.gz" "${sample}/matrix.mtx.gz"
	mv "${sample}_features.tsv.gz" "${sample}/features.tsv.gz"
	mv "${sample}_barcodes.tsv.gz" "${sample}/barcodes.tsv.gz"
done
