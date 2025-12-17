# R scripts

## Introduction

This folder contains the scripts used in R.

## Scripts with descriptions

[accession_filterer.R](accession_filterer.R) - [DEPRECATED] before I found out about DuckDB, I would have to filter the data for quality information the hard way using raw R code, very difficult and messy.

[Create_notebooks.R](Create_notebooks.R) - due to the nature of the way I chose to take notes (i.e. as a quarto book), I wanted to automate some of the process of generating new pages. This works on the command line, using code recorded in the file to create a new page.

[fastafile_getter.R](fastafile_getter.R) - takes a file of accession IDs as an input and runs the NCBI REST API to obtain their sequences, then stores them as [.rds](https://www.sthda.com/english/wiki/saving-data-into-r-data-format-rds-and-rdata) files.

[genomic_fasta_extractor.R](genomic_fasta_extractor.R) - When I get the sequences back from the NCBI, some already have translated proteins, but not all. Thus, it is simplest to just keep the genomic fastas and translate them myself.

[Get_all_pseudo.R](Get_all_pseudo.R) - [DEPRECATED] before DuckDB and using Claude.ai for grouping, when I was still working on just P. aeruginosa, this script basically did the first part of the pipeline in terms of downloading metadata, flattening it, and attempting to group by host.

[GTDB-Tk_tree_maker.R](GTDB-Tk_tree_maker.R) - reads in the newick file output by GTDB-Tk, makes it into a phylogenetic tree with ggtree, annotates with metadata.

[host_metadata_informatics.R](host_metadata_informatics.R) - [DEPRECATED] now done in DuckDB, makes descriptive tables and figures around the host data, not all of them are wholly relevant for the project, it was more data exploration.

[KEGG_enrichment_test_script.R](KEGG_enrichment_test_script.R) - [TEST] script
following the [clusterprofiler book](https://yulab-smu.top/biomedical-knowledge-mining-book/022-kegg.html)
. Testing for how to use enrichKEGG ahead of performing the full analysis.

[metadata_obtaining_script.R](metadata_obtaining_script.R) - [DEPRECATED] old method for downloading metadata with the REST API, complicated and limited to 4,000 output files, too small for current methods. Replaced by manually downloading from the front-end website, as discussed in notebooks (~17-19).

[NCBI_functions.R](NCBI_functions.R) - contains several functions for standardising interactions with the NCBI API, namely attempts at doing API downloads inside of functions. Really, was just an attempt at learning how to do functions in R as it is complecated and interesting, not of great use in modern pipeline.

[pagination_test.R](pagination_test.R) - [TEST] before switch to front-end downloads, tested for using pagination to get over the 4,000 file limit, creates a modified function seen in the previous script.

[reorder_annotations.R](reorder_annotations.R) - [DEPRECATED] before switch to Pseudomonas, was training pipeline on 5 genera, once they exited the annotation stage of the pipeline, they were all in 1 file, needed a script to split into 5 files based on genus.

[supressed_accession_remover.R](supressed_accession_remover.R) - created before switch to Pseudomonas, removes deprecated files from contention as not useful for analysis. Not given deprecated flag as still relevant to do, if there is time.

[test_command_line.R](test_command_line.R) - [TEST] initial test to better understand how the console and terminal worked in Rstudio, and how they could be utilised. In short, the script could be called on the command line, with the argument of rds files. The script would then pull data out of the argument .rds files. This would have been aprototype to the code that went into the [Create_notebooks.R](Create_notebooks.R) script.

[test_fasta_type.R](test_fasta_type.R) - [TEST] early test to understand whether I should use the genomic or protein fasta gotten off of the NCBI. This file would extrat them both, it found that not all of them had protein fastas, so genomic was used

[Utility_functions.R](Utility_functions.R) - commonly used code / useful for good code practice. Script contains 2 functions, 1 for downloading private information needed to access NCBI API, 2 checks if a file path I supply exists, if it does not, the function creates it, avoiding errors.

There are older code scripts in "99_archive" at the top level, but these are even more depricated.
