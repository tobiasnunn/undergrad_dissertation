#!/bin/bash
#SBATCH --account=scw2160
#SBATCH --job-name=TN_5_emapper_search
#SBATCH --partition=compute
#SBATCH --time=02:00:00
#SBATCH --mem=32g
#SBATCH --cpus-per-task=4
#SBATCH -o /scratch/scw2160/09_logs/TN_5_emapper_search_old_params_%j.out
#SBATCH -e /scratch/scw2160/09_logs/TN_5_emapper_search_old_params_%j.err

#increase memory to 32g for the old parameters


# Setup conda environment
eval "$(/apps/languages/anaconda/2024.02/bin/conda shell.bash hook)"
source activate

# Load eggNOG-mapper module
module load eggnog-mapper/2.1.12

# Search phase - new parameters
#time emapper.py \
#  -m diamond \
#  --no_annot \
#  --no_file_comments \
#  -i /scratch/scw2160/TN_area/01_inputs/random_5_combined.faa \
#  --output_dir /scratch/scw2160/TN_area/02_outputs \
#  -o search_5_accessions_new_params \
#  --itype proteins \
#  --sensmode fast \
#  --pident 30 \
#  --query_cover 50 \
#  --evalue 1e-5 \
#  --cpu 4
  
# Search phase - old parameters
time emapper.py \
  -m diamond \
  --no_annot \
  --no_file_comments \
  -i /scratch/scw2160/TN_area/01_inputs/random_5_combined.faa \
  --output_dir /scratch/scw2160/TN_area/02_outputs \
  -o search_5_accessions_old_params \
  --itype proteins \
  --score 60 \
  --pident 40 \
  --query_cover 20 \
  --subject_cover 20 \
  --evalue 0.001 \
  --cpu 4

