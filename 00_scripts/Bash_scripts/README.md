# Bash scripts

## Introduction

This folder contains a few bash scripts used mainly for testing.

## Scripts with descriptions

[emailer.sh](emailer.sh) - this script monitors the `_working` folder where FASTA files get extracted to by [genomic_fasta_extractor.R](../R_scripts/genomic_fasta_extractor.R). When the number of files in the folder changes, it emails me. This lets me monitor the progress of the job since Hawk is no longer sending emails.

[prodigal_test_script.sh](prodigal_test_script.sh) - [TEST] this script was used to test the process of using Prodigal to create protein FASTAs, which is now done through [prodigal_translationiser.R](../Slurm_scripts/prodigal_translationiser.R)

[testing_select_5_fastas.sh](testing_select_5_fastas.sh) - [TEST] this script tested how to combine fasta files together, which is now done through [prodigal_translationiser.R](../Slurm_scripts/prodigal_translationiser.R)