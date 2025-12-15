#!/bin/bash
#SBATCH --account=scw2160
#SBATCH --job-name=annotation_reorderer
#SBATCH --partition=compute           # compute (max resources: 192g; ), highmem (max resources: 384g; 40 threads; 72h)
#SBATCH --time=03:00:00                # Run time
#SBATCH --mem=10g                     # Memory pool for all cores (in MB: 4096 == 4 GB)
#SBATCH -o /scratch/scw2160/09_logs/annotation_reorderer%A.out
#SBATCH -e /scratch/scw2160/09_logs/annotation_reorderer%A.err
#SBATCH --mail-user=0rtpjri2@anonaddy.me
#SBATCH --mail-type=ALL


cd ~/source

module load singularity

directory=/scratch/scw2160/TN_area/02_outputs/

singularity exec ~/tidyverse_latest.sif Rscript 00_scripts/R_scripts/reorder_annotations.R $directory
