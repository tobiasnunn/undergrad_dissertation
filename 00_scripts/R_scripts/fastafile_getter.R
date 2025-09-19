
args <- commandArgs(trailingOnly = TRUE)

#Check if taxon ID was provided
if (length(args) == 0) {
 stop("Please provide a taxonomy ID as an argument")
}

# debug: args <- c(1696)
taxon_ID <- args[1]

# setup hollow for error file
log_file <- paste0("07_log_files/log_genus_", taxon_ID, ".txt")


file_name <- here::here("01_inputs/03_metadata/accession_information_43668_1696_33882_53335_13687_2025-09-10.rds")
all_genera_accession_list <- readRDS(file_name)

accession_list <- all_genera_accession_list[[as.character(taxon_ID)]]



# make sure accessions are not being repeated
file_list <- data.frame(filename = list.files(path = here::here("01_inputs/04_fastas/"), pattern = "*.rds"))
file_list$accession <- gsub(".rds", "", file_list$filename)
accession_list <- accession_list[!names(accession_list) %in% file_list$accession]

log_msg <- paste(Sys.time(), "Number of accessions to process:", length(accession_list))
cat(log_msg, "\n", file = here::here(log_file), append = TRUE)

#debug lines
#random_number <- sample(5:50,3, replace=F) 
#accession_list <- accession_list[random_number]


source(here::here("00_scripts/R_scripts/NCBI_functions.R"))

#run query and return value, save as .rds
#pb <- txtProgressBar(min = 0, max = length(accession_list), style = 3)

for (i in seq_along(accession_list)) {
  tryCatch({
    # i <- 1
    accession <- accession_list[[i]]$accession
    # Args: ftype can be "genome", "protein" or "both"
    fasta_list <- get_dataset_by_accession(accession, "both")
    # debug-line: "GCF_016027095.1"
    saveRDS(fasta_list, file = paste0(here::here("01_inputs/04_fastas/"), accession, ".rds"))
    # log the success
    log_msg <- paste(Sys.time(), "Successfully processed accession",
                       accession) 
    cat(log_msg, "\n", file = here::here(log_file), append = TRUE)
  }, error = function(e) {
    # send error messages to a file
    error_msg <- paste(Sys.time(), "Error with accession",
                       accession, ":", e$message) 
    cat(error_msg, "\n", file = here::here(log_file), append = TRUE)

  })
#  setTxtProgressBar(pb, i)
  Sys.sleep(1)
}
#close(pb)
