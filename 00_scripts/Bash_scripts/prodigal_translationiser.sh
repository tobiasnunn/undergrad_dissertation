#!/bin/bash
#SBATCH --account=scw2160
#SBATCH --job-name=convert_with_prodigal
#SBATCH --partition=compute           # compute (max resources: 192g; ), highmem (max resources: 384g; 40 threads; 72h)
#SBATCH --time=4:00:00                # Run time
#SBATCH --mem=8g                     # Memory pool for all cores (in MB: 4096 == 4 GB)
#SBATCH -o /scratch/scw2160/09_logs/prodigal_translationiser_%A.out
#SBATCH -e /scratch/scw2160/09_logs/prodigal_translationiser_%A.err

# Load the Prodigal module
module load Prodigal/2.6.3

#path to where .fna files are
cd ~/source/01_inputs/04_fastas/_working

# Create list of accessions that need processing
echo "Checking for files that need processing..."
to_process=()
for file in *.fna; do
    if [[ -f "$file" ]]; then
        full_name=$(basename "$file" .fna)
        accession=$(echo "$full_name" | sed 's/_genome$//')
        
        # Check if already processed
        if [[ ! -f "${accession}_renamed.faa" ]]; then
            to_process+=("$file")
            echo "Will process: $file -> $accession"
        else
            echo "Already processed: $file (${accession}_renamed.faa exists)"
        fi
    fi
done

echo "Found ${#to_process[@]} files to process"

# Process files
for file in "${to_process[@]}"; do
    full_name=$(basename "$file" .fna)
    accession=$(echo "$full_name" | sed 's/_genome$//')
    
    echo "Processing $file, accession: $accession"
    
    # Run Prodigal
    prodigal -i "$file" -a "${accession}_proteins.faa" > /dev/null 2>&1
    
    # Replace contig identifiers
    awk -v acc="$accession" '/^>/ {gsub(/^>[^_]*_[0-9]*/, ">"acc"_"++counter)} 1' "${accession}_proteins.faa" > "${accession}_renamed.faa"
    
    # Clean up intermediate files
    rm "$file" "${accession}_proteins.faa"
    echo "Cleaned up $file and ${accession}_proteins.faa"
done

# Combine all renamed files (including any that existed before)
echo "Combining all *_renamed.faa files..."
cat *_renamed.faa > combined_pseudomonas_proteins.faa
echo "Created combined_proteins.faa with $(grep -c "^>" combined_pseudomonas_proteins.faa) total proteins"