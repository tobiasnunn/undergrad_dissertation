# list rds files and extract genome fastas to make ready for eggnog search

source(here::here("00_scripts/R_scripts/Utility_functions.R"))
check_and_create_dir("01_inputs/04_fastas/_working")

# get list of rds files
file_name <- here::here("01_inputs/03_metadata/filtered_accession.tsv")
all_genera_accession_list <- read.delim(file_name)
all_genera_accession_list$filepath <- paste0("01_inputs/04_fastas/", all_genera_accession_list$accession, ".rds")
print(all_genera_accession_list$filepath)

# read list in, in order to extract genome
# extract genome file from each read in accession, update contig name
#TODO: make sure this does not pad with NAs if the list is too short
#process_list <- file.path("01_inputs/04_fastas", paste0(total_fasta_list[201:length(total_fasta_list)]))
process_list <- all_genera_accession_list$filepath

for (rds_file in process_list) {
  tryCatch({
      # rds_file <- process_list[1]
      # Load the RDS file
      fasta_data <- readRDS(rds_file)
      
      # Extract accession from filename
      accession <- tools::file_path_sans_ext(basename(rds_file))
      
      # Extract genome fasta
      genome_fasta <- fasta_data[[accession]]$genome  
      
      # Update header lines to include accession at the beginning
      #header_lines <- grepl("^>", protein_fasta)
      #protein_fasta[header_lines] <- paste0(">", accession, "_", substring(protein_fasta[header_lines], 2))
      
      # Save genome fasta to file
      output_file <- file.path("01_inputs/04_fastas/_working", paste0(accession, "_genome.fna"))
      writeLines(genome_fasta, output_file)  
      print(paste0("Successfully extracted the genome for ", accession, 
                   ". File path = ", output_file))
    }, error = function(e) {
      # print error messages
      error_msg <- paste(Sys.time(), "Error with accession",
                       accession, ":", e$message) 
      print(error_msg)

    })
}
print(Sys.time())
sink(file = NULL)