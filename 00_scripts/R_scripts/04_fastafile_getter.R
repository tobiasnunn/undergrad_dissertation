
file_name <- "../01_inputs/03_metadata/accession_information_43668_1696_33882_53335_13687_2025-09-10.rds"
accession_list <- readRDS(file_name)

source("../00_scripts/R_scripts/NCBI_functions.R")

#testing for loop
for (i in 1:5) {
  # Args: ftype can be "genome", "protein" or "both"
  fasta_list <- get_dataset_by_accession(accession_list$`1696`[[i]]$current_accession, "both")
  # debug-line: "GCF_016027095.1"
  saveRDS(fasta_list, file = paste0("../01_inputs/04_fastas/", accession_list$`1696`[[i]]$current_accession, ".rds"))
}
