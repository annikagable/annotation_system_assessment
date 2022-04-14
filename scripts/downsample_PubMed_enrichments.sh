#!/bin/bash

## This script will copy all enrichment results to a new folder, 
## and in the process downsample the PubMed results so that the 
## number of results is comparable to other annotation systems
## and easier to plot for the term size vs effect size plot.

data_dir="$1"
sampling_factor="$2"
out_dir="$3"

mkdir -p "$out_dir"

for f in "$data_dir"/*; do
    out_f="$out_dir"/$(basename "$f")
    if [[ "$f" =~ "PubMed" ]]
    then
        n_lines=$(< "$f" wc -l)
        head -n1 "$f" > "$out_f"
        tail -n+2 "$f" | shuf -n $((n_lines/sampling_factor)) >> "$out_f" 
    else
        cp "$f" "$out_f"
    fi
done
