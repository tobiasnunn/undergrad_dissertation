#!/bin/bash
#SBATCH --account=scw2160
#SBATCH --job-name=TN_1000_annotation
#SBATCH --partition=highmem
#SBATCH --time=06:00:00
#SBATCH --mem=50g
#SBATCH --cpus-per-task=4
#SBATCH -o /scratch/scw2160/09_logs/TN_1000_annotation_%j.out
#SBATCH -e /scratch/scw2160/09_logs/TN_1000_annotation_%j.err

source activate
module load eggnog-mapper/2.1.12

time emapper.py \
  --annotate_hits_table /scratch/scw2160/TN_area/02_outputs/search_1000_accessions.emapper.seed_orthologs \
  --no_file_comments \
  --output_dir /scratch/scw2160/TN_area/02_outputs \
  -o annotation_1000_accessions \
  --dbmem \
  --cpu 4
