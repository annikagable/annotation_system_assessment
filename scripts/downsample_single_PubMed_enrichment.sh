#!/bin/bash

## This script will copy an enrichment result for a particular user
## input to a new folder, and in the process downsample the PubMed 
## results so that the number of results is comparable to other 
## annotation systems and thus easier to plot for the 
## 'term size vs effect size' plot.

in_f="$1"
sampling_factor="$2"
out_f="$3"

mkdir -p $(dirname "$out_f")

if [[ "$in_f" =~ "PubMed" ]]
then
    n_lines=$(< "$in_f" wc -l)
    head -n1 "$in_f" > "$out_f"
    tail -n+2 "$in_f" | shuf -n $((n_lines/sampling_factor)) >> "$out_f" 
else
    cp "$in_f" "$out_f"
fi
