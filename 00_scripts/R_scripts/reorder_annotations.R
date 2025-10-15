# out of EggNOG the annotations are in abstract batches
# This script separates them out by genus after combining them

library(tidyverse)

# 1. read in metadata file ------------------------------------------------
# so that I have a list to compare against later

accession_info_by_genus <- readRDS(here::here("01_inputs/03_metadata/accession_information_43668_1696_33882_53335_13687_2025-09-10.rds"))

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

rm(accession_info_by_genus)
# 2. read in gzip files ---------------------------------------------------

KO_files <- list.files("C:/Users/tobyn/OneDrive/Work/tnunn_research/2025 file storage/", 
                       pattern = "*.annotations.copy.gz", full.names = T)

KO_file_list <- readr::read_tsv(KO_files, id = "file_name")


KO_file$accession <- gsub("_[0-9]+$", "", KO_file$`#query`)
first_sample <- filter(KO_file, accession == "GCA_016429035.1")
unique(KO_file$accession)


# 3. code for hawk --------------------------------------------------------

accession_info_by_genus <- readRDS(here::here("01_inputs/03_metadata/accession_information_43668_1696_33882_53335_13687_2025-09-10.rds"))

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

