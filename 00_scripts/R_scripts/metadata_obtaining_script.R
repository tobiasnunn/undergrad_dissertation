
# 01. libraries and loading custom functions -----------------------------------------------------------

# path relative to the .proj file
source(here::here("00_scripts/R_scripts/NCBI_functions.R"))


# 02. obtain taxon IDs for the 5 genera -----------------------------------

# call function to obtain tax IDs from NCBI

# vector of genus names
genera <- c("Brachybacterium", "Brevibacterium", "Microbacterium", "Pantoea", "Sphingomonas")

taxname_with_ID <- get_taxID(genera, "genus")

# 03. download metadata (using cool functions) ---------------------------------------------------


result <- tryCatch({
  download_taxonomic_information_by_tax_ID(taxname_with_ID$tax_id, here::here("01_inputs/03_metadata/"))
}, error = function(e) {
  message("Error in download_metadata: ", e$message)
  stop("Download failed: ", e$message)
})

# check to make sure they were written right
banana <- readRDS(here::here("01_inputs/03_metadata/taxonomy_information_Sphingomonas_Brevibacterium_Microbacterium_Brachybacterium_Pantoea_2025-09-10.rds"))
coconut <- readRDS(here::here("01_inputs/03_metadata/accession_information_43668_1696_33882_53335_13687_2025-09-10.rds"))
