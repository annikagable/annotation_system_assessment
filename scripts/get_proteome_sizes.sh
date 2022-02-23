#!/bin/bash
# This script takes all stringIds from STRING's protein.info file, splits the taxIs off, and counts the number
# of lines with the same taxId. The output are two columns: 1) the taxId and 2) number of protein-coding genes.
#
# Usage: get_proteome_sizes.sh protein.info.v11.0.txt.gz > taxId_to_proteome_size.tsv 
   

content=$(zcat $1 | tail -n+2 | cut -f1 | cut -d'.' -f1 | uniq -c | tr -s ' ' '\t'| awk 'OFS="\t" {print $2, $1}')
header="taxId\tproteome_size"
echo -e "$header\n$content"
