#!/bin/bash
#SBATCH --account=scw2160
#SBATCH --job-name=convert_with_prodigal
#SBATCH --partition=compute           # compute (max resources: 192g; ), highmem (max resources: 384g; 40 threads; 72h)
#SBATCH --time=10:00:00                # Run time
#SBATCH --mem=8g                     # Memory pool for all cores (in MB: 4096 == 4 GB)
#SBATCH -o /scratch/scw2160/09_logs/prodigal_translationiser_%A.out
#SBATCH -e /scratch/scw2160/09_logs/prodigal_translationiser_%A.err
#SBATCH --mail-user=0rtpjri2@anonaddy.me
#SBATCH --mail-type=ALL

echo "Job $SLURM_JOB_ID started at $(date)" | \
    mail -s "Job Started: $SLURM_JOB_ID" 0rtpjri2@anonaddy.me

# Load the Prodigal module
module load Prodigal/2.6.3

#path to where .fna files are
cd ~/source/01_inputs/04_fastas/_running

# Run the main processing in background
{
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

  cat *_renamed.faa > 1000_1_combined_proteins.faa
} &
main_job_pid=$!

# Monitor progress in the background
{
    while kill -0 $main_job_pid 2>/dev/null; do
        sleep 1200  # 20 minutes
        if kill -0 $main_job_pid 2>/dev/null; then
            faa_count=$(ls -1 *_renamed.faa 2>/dev/null | wc -l)
            echo "Job $SLURM_JOB_ID: $faa_count .faa files processed at $(date)" | \
                mail -s "Prodigal Progress: $faa_count files" 0rtpjri2@anonaddy.me
        fi
    done
} &

# Wait for main job to finish
wait $main_job_pid

# Final notification
final_count=$(ls -1 *_renamed.faa 2>/dev/null | wc -l)
echo "Job $SLURM_JOB_ID completed with $final_count .faa files at $(date)" | \
    mail -s "Prodigal Complete: $final_count files" 0rtpjri2@anonaddy.me
