#! /bin/bash

## From the protein.info.txt.gz file, extract each species into a file of its own.
## Example ./split_protein_info_by_species.sh  protein.info.v11.0.txt.gz out_dir

zipped=$1
dir=$2
unzipped_base=$(basename $zipped .gz)

zcat $zipped | tail -n+2 | awk -v dir=$dir -v base=$unzipped_base 'BEGIN {FS="."} {print > dir"/"$1"."base}'

# Explainer: 
# unzip; remove header; pass two variables to awk; use dot as field separator; print 
# content of file to a file whose name consists of the first field, 
# i.e. the species, and base name of original file. 
