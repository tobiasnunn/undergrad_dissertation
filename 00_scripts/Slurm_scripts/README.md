# Slurm scripts

## Introduction
This folder contains the scripts used by [Slurm](https://slurm.schedmd.com/documentation.html) on [SuperComputing Wales](https://www.supercomputing.wales/) high performance computing environment.

## Scripts with descriptions
[eggnog_annotator.sh](eggnog_annotator.sh) - starts the second phase of the eggNOG-mapper annotation process. It runs on the high mem partition on Hawk because the entire eggNOG database (c. 44GB) has to be loaded into memory.

[eggnog_searcher.sh](eggnog_searcher.sh) - performs the first phase of the eggNOG-mapper annotation process, performing a DIAMOND search against the combined protein FASTA created by [prodigal_translationiser.sh](prodigal_translationiser.sh).

[fasta_getter_runner.sh](fasta_getter_runner.sh) - uses Slurm to execute the R script [fastafile_getter.R](../R_scripts/fastafile_getter.R).

[gtdbtk_runner.sh](gtdbtk_runner.sh) - runs the [gtdb-tk de novo workflow](https://ecogenomics.github.io/GTDBTk/commands/de_novo_wf.html) on a set of genomic fastas downloaded from the NCBI via [fasta_getter_runner.sh](fasta_getter_runner.sh). It stores the GTDB-TK output in /source/02_analysis/02_gtdbtk.

[prodigal_translationiser.sh](prodigal_translationiser.sh) - looks for a set of genomic FASTAs `(*.fna)` in 01_inputs/04_fastas/_working and uses [Prodigal](https://github.com/hyattpd/Prodigal) 2.6.3 to create protein FASTAs `(*.faa)`. Because the protein FASTAs all have to get combined into one big file for the eggNOG-mapper processes, the script also updates the contig labels to include the accession ID. This ensures the eggNOG annotations can be attributed to the right accession at the end.

[reorder_annotations_runner.sh](reorder_annotations_runner.sh) - uses Slurm to execute the R script [reorder_annotations.R](../R_scripts/reorder_annotations.R).

