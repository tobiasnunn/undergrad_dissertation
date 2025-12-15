# file for removing supressed accessions (should become redundant fairly quick)

file_name <- here::here("01_inputs/03_metadata/accession_information_43668_1696_33882_53335_13687_2025-09-10.rds")
all_genera_accession_list <- readRDS(file_name)

# obtain list of supressed accessions
suppressed_accessions_list <- c()  

for (genus_id in names(all_genera_accession_list)) {
  genus_data <- all_genera_accession_list[[genus_id]]
  for (acc_id in names(genus_data)) {
    if (genus_data[[acc_id]]$assembly_info$assembly_status == "suppressed") {
      suppressed_accessions_list <- c(suppressed_accessions_list, acc_id)
    }
  }
}

# delete supressed accessions from all_genera_accession_list
## Construct the full file paths
fasta_files_to_delete <- file.path("01_inputs/04_fastas", paste0(suppressed_accessions_list, ".rds"))

## Check which files actually exist before trying to delete
existing_files <- fasta_files_to_delete[file.exists(fasta_files_to_delete)]

## Delete the existing files
if (length(existing_files) > 0) {
  file.remove(existing_files)
  cat("Deleted", length(existing_files), "suppressed accessions\n")
} else {
  cat("No suppressed accessions found to delete\n")
}