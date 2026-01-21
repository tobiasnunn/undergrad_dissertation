# file for reading in list of accessions made elsewhere and downloading fasta files

file_name <- "01_inputs/03_metadata/filtered_accession.tsv"
accession_list <- read.delim(file_name)

# make sure accessions are not being repeated
file_list <- data.frame(filename = list.files(path = "01_inputs/04_fastas/"), pattern = "*.rds")
file_list$accession <- gsub(".rds", "", file_list$filename)
accession_list <- accession_list[!accession_list %in% file_list$accession]

log_msg <- paste(Sys.time(), "Number of accessions to process:", length(accession_list$accession))
print(log_msg)

source("00_scripts/R_scripts/NCBI_functions.R")

#run query and return value, save as .rds
#pb <- txtProgressBar(min = 0, max = length(accession_list), style = 3)

for (i in seq_along(accession_list$accession)) {
  tryCatch({
    # i <- 1
    #accession <- accession_list[[i]]$accession
    accession <- accession_list$accession[i]
    # Args: ftype can be "genome", "protein" or "both"
    fasta_list <- get_dataset_by_accession(accession, "genome")
    # debug-line: "GCF_016027095.1"
    saveRDS(fasta_list, file = paste0("01_inputs/04_fastas/"), accession, ".rds")
    # log the success
    log_msg <- paste(Sys.time(), "Successfully processed accession",
                       accession) 
    print(log_msg)
  }, error = function(e) {
    # send error messages to a file
    error_msg <- paste(Sys.time(), "Error with accession",
                       accession, ":", e$message) 
    print(error_msg)

  })
#  setTxtProgressBar(pb, i)
  Sys.sleep(1)
}
#close(pb)
