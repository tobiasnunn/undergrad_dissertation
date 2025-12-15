# out of EggNOG the annotations are in abstract batches
# This script separates them out by genus after combining them

library(tidyverse)

# # 1. read in metadata file ------------------------------------------------
# # so that I have a list to compare against later
# 
# accession_info_by_genus <- readRDS(here::here("01_inputs/03_metadata/accession_information_43668_1696_33882_53335_13687_2025-09-10.rds"))
# 
# ## MK 4 - tidyverse, without base R (Claude assisted)
# # code adapted from "host_metadata_informatics.R"
# genus_with_accession <- imap(accession_info_by_genus, function(taxon_data, tax_id){
#   imap(taxon_data, function(accession_data, accession_id){
#     data.frame(
#       tax_id = as.integer(tax_id),
#       accession_id = accession_id
#      # host = pluck(accession_data, "assembly_info", "biosample", "host", .default = NA),
#      # isolation_source = pluck(accession_data, "assembly_info", "biosample", "isolation_source", .default = NA)
#     )
#   }) %>% list_rbind()
# }) %>% list_rbind()
# 
# rm(accession_info_by_genus)
# # 2. read in gzip files ---------------------------------------------------
# 
# KO_files <- list.files("C:/Users/tobyn/OneDrive/Work/tnunn_research/2025 file storage/", 
#                        pattern = "*.annotations.copy.gz", full.names = T)
# 
# KO_file_list <- readr::read_tsv(KO_files, id = "file_name")
# 
# 
# KO_file$accession <- gsub("_[0-9]+$", "", KO_file$`#query`)
# first_sample <- filter(KO_file, accession == "GCA_016429035.1")
# unique(KO_file$accession)
# 
# # test split()
# bazinga <- data.frame(name = c("banana", "bat", "bannana!", "cat", "joe stinky butt", "coffee"), 
#                       category = c(69, 1, 356, 69, 69, 1))
# bazinga2 <- split(bazinga, bazinga$category)

# 3. code for hawk --------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

#Check if taxon ID was provided
if (length(args) == 0) {
  stop("Please provide the directory the zip files are located in.")
}

# debug: directory <- here::here("01_inputs/05_annotations/original")
directory <- args[1]

# load utility function so script can create output site
source(here::here("00_scripts/R_scripts/Utility_functions.R"))
check_and_create_dir(here::here("01_inputs/05_annotations/by_genus"))
  
# setup hollow for error file
log_file <- paste0("07_log_files/log_annotation_", Sys.Date(), ".txt")


# read in total list of accessions by genus for comparison
accession_info_by_genus <- readRDS(here::here("01_inputs/03_metadata/accession_information_43668_1696_33882_53335_13687_2025-09-10.rds"))

# make it into a dataframe so I can compare against it
## MK 4 - tidyverse, without base R (Claude assisted)
# code adapted from "host_metadata_informatics.R"
genus_with_accession <- imap(accession_info_by_genus, function(taxon_data, tax_id){
  imap(taxon_data, function(accession_data, accession_id){
    data.frame(
      tax_id = as.integer(tax_id),
      accession_id = accession_id
      # host = pluck(accession_data, "assembly_info", "biosample", "host", .default = NA),
      # isolation_source = pluck(accession_data, "assembly_info", "biosample", "isolation_source", .default = NA)
    )
  }) %>% list_rbind()
}) %>% list_rbind()

# read in list of files so the for loop has smth to iterate over
KO_files <- list.files(here::here(directory), 
                       pattern = "*.annotations.copy.gz", full.names = T)
print(paste0("Found ", length(KO_files), " zip files in the supplied directory."))
print(paste0("KO files:", KO_files, collapse = ","))

# for loop to read in each file in turn, lookup the genus in "genus_with_accession" and copy the relevant lines
# to a new dataframe specific to each genus, so the thing will have essentially been split by genus
for (KO_file in KO_files) {
  #debug: KO_file <- KO_files[1]
  # read in file
  annotations <- readr::read_tsv(KO_file, id = "file_name") %>% 
  #rename important column with weird name
    rename(query = `#query`)
  
  # make new column with just accession (cut the contig part of the "query")
  annotations$accession <- substr(annotations$query, 1, 15)
  
  # combine the two objects so genus information transfers to "annotations"
  annotations <- annotations %>% 
    left_join(genus_with_accession, by = join_by(accession == accession_id))
  
  # chop up annotations by tax_id
  split_annotations <- split(annotations, annotations$tax_id)
  # debug
  print(paste0("number of genera found : ", length(split_annotations)))
  
  # store as individual RDS files
  for (j in seq_along(split_annotations)) {
    # read in genus code names so that RDS files have good names
    rds_name <- paste0("01_inputs/05_annotations/by_genus/annotations_",names(split_annotations)[j], ".rds")
    
    # grab object in position j
    annotation_df <- split_annotations[[j]]
    # debug
    print(paste0("number of rows in ", rds_name, ": ", nrow(annotation_df)))
    
    # check to see if the outer loop is on first iteration or not
    if(file.exists(here::here(rds_name))){
      #monitoring
      print(paste0(rds_name, " already exists. Adding to file..."))
      # if the thing exists, I want to read it in, to combine it with new material
      existing_df <- readRDS(here::here(rds_name))
      #combine existing_df and annotation_df
      annotation_df <- bind_rows(existing_df, annotation_df)
      }
    
    #save annotation_df out as rds with rds_name as file name
    saveRDS(annotation_df, file = paste0(here::here(rds_name)))
  }
}
