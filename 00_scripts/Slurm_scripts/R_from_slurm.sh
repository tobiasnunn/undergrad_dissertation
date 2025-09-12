#!/bin/bash
#SBATCH --account=scw2160
#SBATCH --job-name=R_test
#SBATCH --partition=compute           # compute (max resources: 192g; ), highmem (max resources: 384g; 40 threads; 72h)
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --time=00:10:00                # Run time
#SBATCH --mem=1g                     # Memory pool for all cores (in MB: 4096 == 4 GB)
#SBATCH -o /scratch/scw2160/09_logs/R_test%A.out
#SBATCH -e /scratch/scw2160/09_logs/R_test%A.err
#SBATCH --mail-user=0rtpjri2@anonaddy.me
#SBATCH --mail-type=ALL

module load R/4.4.2

Rscript ~/source/00_scripts/R_scripts/test_command_line.R ~/source/01_inputs/03_metadata/taxonomy_information_Sphingomonas_Brevibacterium_Microbacterium_Brachybacterium_Pantoea_2025-09-10.rds
