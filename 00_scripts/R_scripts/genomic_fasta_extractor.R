# list rds files and extract protein fastas to make ready for eggnog search
#NOTE: this is an innitial build with some stuff (resumability mainly) left to add, 
# just for the test 200 to begin with

source(here::here("00_scripts/R_scripts/Utility_functions.R"))
check_and_create_dir("01_inputs/04_fastas/_working")

# get list of rds files
if(!file.exists(here::here("01_inputs/04_fastas/total_fasta_list.txt"))){
  file_list <- list.files(path = "01_inputs/04_fastas/", pattern = "*.rds")
  writeLines(file_list, here::here("01_inputs/04_fastas/total_fasta_list.txt"))
}

# read list in, in order to extract protein
total_fasta_list <- readLines(here::here("01_inputs/04_fastas/total_fasta_list.txt"))

# extract protein file from each read in accession, update contig name
#TODO: make sure this does not pad with NAs if the list is too short
process_list <- file.path("01_inputs/04_fastas", paste0(total_fasta_list[1:200]))

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
      
      # Save protein fasta to file
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
