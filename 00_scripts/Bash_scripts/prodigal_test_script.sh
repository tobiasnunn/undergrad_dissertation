#!/bin/bash

# Load the Prodigal module
module load Prodigal/2.6.3

#path to where .fna files are
cd ~/source/01_inputs/04_fastas/_working

# Run prodigal on all .fna files
for file in *.fna; do
    full_name=$(basename "$file" .fna)
    # Extract just the accession (remove _genome suffix)
    accession=$(echo "$full_name" | sed 's/_genome$//')
    echo "Processing $file, accession: $accession"
    prodigal -i "$file" -a "${accession}_proteins.faa"
    # Replace the contig identifier with the accession number
    awk -v acc="$accession" '/^>/ {gsub(/^>[^_]*_[0-9]*/, ">"acc"_"++counter)} 1' "${accession}_proteins.faa" > "${accession}_renamed.faa"
done

cat *_renamed.faa > combined_proteins.faa
