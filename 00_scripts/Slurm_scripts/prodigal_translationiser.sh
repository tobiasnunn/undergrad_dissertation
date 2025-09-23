#!/bin/bash
#SBATCH --account=scw2160
#SBATCH --job-name=convert_with_prodigal
#SBATCH --partition=compute           # compute (max resources: 192g; ), highmem (max resources: 384g; 40 threads; 72h)
#SBATCH --time=3:00:00                # Run time
#SBATCH --mem=8g                     # Memory pool for all cores (in MB: 4096 == 4 GB)
#SBATCH -o /scratch/scw2160/09_logs/prodigal_translationiser_%A.out
#SBATCH -e /scratch/scw2160/09_logs/prodigal_translationiser_%A.err

# Load the Prodigal module
module load Prodigal/2.6.3

#path to where .fna files are
cd ~/source/01_inputs/04_fastas/_running_1500

# Run prodigal on all .fna files
for file in *.fna; do
    full_name=$(basename "$file" .fna)
    # Extract just the accession (remove _genome suffix)
    accession=$(echo "$full_name" | sed 's/_genome$//')
    echo "Processing $file, accession: $accession"
    prodigal -i "$file" -a "${accession}_proteins.faa" > /dev/null 2>&1
    # Replace the contig identifier with the accession number
    awk -v acc="$accession" '/^>/ {gsub(/^>[^_]*_[0-9]*/, ">"acc"_"++counter)} 1' "${accession}_proteins.faa" > "${accession}_renamed.faa"
done

cat *_renamed.faa > 1500_1_combined_proteins.faa