#!/bin/bash
#SBATCH --account=scw2160
#SBATCH --job-name=fasta_downloader
#SBATCH --partition=compute           # compute (max resources: 192g; ), highmem (max resources: 384g; 40 threads; 72h)
#SBATCH --time=02:00:00                # Run time
#SBATCH --mem=1g                     # Memory pool for all cores (in MB: 4096 == 4 GB)
#SBATCH -o /scratch/scw2160/09_logs/fasta_downloader%A.out
#SBATCH -e /scratch/scw2160/09_logs/fasta_downloader%A.err
#SBATCH --mail-user=0rtpjri2@anonaddy.me
#SBATCH --mail-type=ALL
#SBATCH --array=1-5

cd ~/source

module load singularity

genus_ids=(43668 1696 13867 33882 53335)
genus_id=${genus_ids[$SLURM_ARRAY_TASK_ID-1]}

singularity exec ~/tidyverse_latest.sif Rscript 00_scripts/R_scripts/fastafile_getter.R $genus_id
