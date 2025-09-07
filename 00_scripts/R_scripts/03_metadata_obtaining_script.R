
# 01. libraries and loading custom functions -----------------------------------------------------------

# path relative to the .proj file
source("../00_scripts/R_scripts/NCBI_functions.R")


# 02. obtain taxon IDs for the 5 genera -----------------------------------

# call function to obtain tax IDs from NCBI

# vector of genus names
genera <- c("Brachybacterium", "Brevibacterium", "Microbacterium", "Pantoea", "Sphingomonas")

taxname_with_ID <- get_taxID(genera, "genus")

# 03. download metadata ---------------------------------------------------

result <- tryCatch({
  download_NCBI_metadata_by_taxID(taxname_with_ID$tax_id, "../01_inputs/03_metadata/")
}, error = function(e) {
  message("Error in download_metadata: ", e$message)
  stop("Download failed: ", e$message)
})
