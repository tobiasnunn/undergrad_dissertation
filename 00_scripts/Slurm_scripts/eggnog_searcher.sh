#!/bin/bash
#SBATCH --account=scw2160
#SBATCH --job-name=TN_2040_emapper_search
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --mem=64g
#SBATCH --cpus-per-task=4
#SBATCH -o /scratch/scw2160/09_logs/TN_2040_emapper_search_%j.out
#SBATCH -e /scratch/scw2160/09_logs/TN_2040_emapper_search_%j.err


# Setup conda environment
eval "$(/apps/languages/anaconda/2024.02/bin/conda shell.bash hook)"
source activate

# Load eggNOG-mapper module
module load eggnog-mapper/2.1.12

# Search phase - new parameters
time emapper.py \
  -m diamond \
  --no_annot \
  --no_file_comments \
  -i /scratch/scw2160/TN_area/01_inputs/combined_2040_proteins.faa \
  --output_dir /scratch/scw2160/TN_area/02_outputs \
  -o search_2040_accessions \
  --itype proteins \
  --sensmode fast \
  --pident 30 \
  --query_cover 50 \
  --evalue 1e-5 \
  --cpu 4
  

