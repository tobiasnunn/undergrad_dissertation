# these are functions for when I interface with the NCBI database


# 01. get taxonomy ID -----------------------------------------------------
# e.g. I want to work on the genus "Sphingomonas", however, I cannot use that name
# to find things in the database, thus, I must find its accompany-ing ID
# I could manually look this up, but that is not clever

get_taxID <- function(taxonomic_name, desired_rank = "genus"){
  # throw an error if "taxonomic_name" is not a character
  stopifnot(
    "taxonomic_name must be character" = is.character(taxonomic_name),
    "taxonomic_name cannot be empty" = length(taxonomic_name) > 0,
    "desired_rank not in list of known ranks" = desired_rank %in% 
      c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"))
  
  source("../00_scripts/R_scripts/Utility_functions.R")
  ncbi_secrets <- get_secrets("ncbi")
  
  genus_info <- data.frame(tax_id = numeric(), sci_name = character()) # setup hollow df
  
  for (taxon in taxonomic_name) {
    #debug
    # taxon <- "Mustelidae" / desired_rank <- "family" / taxon <- "Homo sapiens"
    message("Processing:", taxon, "\n")
    
    encoded_taxon <- utils::URLencode(taxon, reserved = TRUE) # error handling for spaces
    
    # connect to NCBI database, bring down data into "req" object
    req <- httr2::request(ncbi_secrets$host) |> 
      httr2::req_auth_bearer_token(ncbi_secrets$api_key) |> 
      httr2::req_url_path_append('taxonomy', 'taxon_suggest', encoded_taxon) |> 
      httr2::req_headers(Accept = 'application/json') |> 
      httr2::req_url_query(tax_rank_filter = "higher_taxon") |> 
      httr2::req_perform()
    
    # store result
    if(httr2::resp_is_error(req)){
      stop(paste("API call failed with error number: ", resp_status()))
    } else{
      result <- httr2::resp_body_json(req)$sci_name_and_ids |>
        purrr::keep(~ !is.null(.x$rank) && .x$rank == toupper(desired_rank)) |> # using the purrr package to filter the list
        purrr::map_dfr(~ data.frame(tax_id = as.numeric(.x$tax_id), sci_name = .x$sci_name))
      
      genus_info <- rbind(genus_info, result) # bind the new row into the set_up dataframe
    }
  }
  return(genus_info)
} 
# 02. get metadata by taxon ID --------------------------------------------
# e.g. I run the first function to get the tax ID, then I can use that
# to lookup all accessions related to that ID (i.e. all Sphingomonas accessions)

# another NCBI REST API, with different params in the "url_path_append"

download_NCBI_metadata_by_taxID <- function(tax_ID, file_location){
  # throw an error if "taxonomic_name" is not a character
  stopifnot(
    "taxonomic_ID must be integer" = is.numeric(tax_ID),
    "taxonomic_ID cannot be empty" = length(tax_ID) > 0,
    "file_location cannot be empty" = length(file_location) > 0)
  
  library(tidyverse)
  
  # check directory "file_location" exists, if not, create
  source("../00_scripts/R_scripts/Utility_functions.R")
  check_and_create_dir(file_location)
  ncbi_secrets <- get_secrets("ncbi")
  
  
  # get the name from the ID for naming outputs better
  #tax_ID <- c(1696, 13687, 33882)
  taxon_info <- data.frame(tax_ID = tax_ID)
  tax_ID_string <- paste0(taxon_info$tax_ID, collapse = ",")
  
  req <- httr2::request(ncbi_secrets$host) |> 
    httr2::req_auth_bearer_token(ncbi_secrets$api_key) |> 
    httr2::req_url_path_append('taxonomy', 'taxon', tax_ID_string) |> 
    httr2::req_headers(Accept = 'application/json') |> 
    httr2::req_url_query(page_size = 1000) |> 
    httr2::req_perform()
  
  # store result
  if(httr2::resp_is_error(req)){
    stop(paste("could not find names for IDs. Error number: ", resp_status()))
  } else{
    resps <- httr2::resp_body_json(req)
    taxonomy_information <- resps$taxonomy_nodes %>%
      map_dfr(~ {
        # Skip nodes with errors
        if (!is.null(.x$errors)) {
          return(NULL)
        }
        
        # Skip if no taxonomy data
        if (is.null(.x$taxonomy)) {
          return(NULL)
        }
        
        tax_data <- .x$taxonomy
        # tax_data <- taxonomy_information[[1]]$taxonomy
        
        # Extract and pivot counts to wide format
        counts_info <- tax_data$counts %>%
          map_dfr(~ tibble(type = .x$type, count = .x$count)) %>%
          pivot_wider(names_from = type, values_from = count, values_fill = 0)
        
        # Combine with other taxonomy fields
        tibble(
          tax_id = tax_data$tax_id %||% NA,
          organism_name = tax_data$organism_name %||% NA,
          blast_name = tax_data$blast_name %||% NA,
          rank = tax_data$rank %||% NA,
          genomic_moltype = tax_data$genomic_moltype %||% NA,
          lineage = list(tax_data$lineage),  # Keep as list column
          children = list(tax_data$children) # Keep as list column
        ) %>%
          bind_cols(counts_info)  # Add the flattened count columns
      })
    saveRDS(taxonomy_information, file = paste0(file_location, "taxonomy_information.rds"))
  }
  # file_location <- "../01_inputs/03_metadata/"
  
  for (taxon in tax_ID) {
    # taxon <- 43668 / file_path <- "../01_inputs/03_metadata_copy/"
    req <- httr2::request(ncbi_secrets$host) |>
      httr2::req_auth_bearer_token(ncbi_secrets$api_key) |>
      httr2::req_url_path_append('genome', 'taxon', taxon, 'dataset_report') |>
      httr2::req_url_query(
        page_size = 200,
        filters.exclude_paired_reports = TRUE) |>
      httr2::req_headers(Accept = 'application/json')

    # keeps calling req until all pages are returned ("page_size  200")
    resps <- httr2::req_perform_iterative(
      req,
      next_req = httr2::iterate_with_cursor(
        'page_token',
        resp_param = function(resp) {
          httr2::resp_body_json(resp)$next_page_token
        }
      ),
      path =  paste0(file_location, taxon, "_{Sys.Date()}_page{i}.json") 
    )
    
    message(paste(taxon, " has been processed"))
  }
  message("all done successfully")
}
