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

download_taxonomic_information_by_tax_ID <- function(tax_ID, file_location){
  # throw an error if "taxonomic_name" is not a character
  stopifnot(
    "taxonomic_ID must be integer" = is.numeric(tax_ID),
    "taxonomic_ID cannot be empty" = length(tax_ID) > 0,
    "file_location cannot be empty" = length(file_location) > 0)
  
  # check directory "file_location" exists, if not, create
  source("../00_scripts/R_scripts/Utility_functions.R")
  check_and_create_dir(file_location) # create dir if not already exist
  ncbi_secrets <- get_secrets("ncbi") # get secrets for API call
  
  # convert vector of tax IDs into single string
  tax_ID_string <- paste0(tax_ID, collapse = ",")
  # tax_ID_string <- c(1696, 33882) |> paste0(collapse = ",")
  
  # query NCBI database using REST API
  req <- httr2::request(ncbi_secrets$host) |> 
    httr2::req_auth_bearer_token(ncbi_secrets$api_key) |> 
    httr2::req_url_path_append('taxonomy', 'taxon', tax_ID_string) |> 
    httr2::req_headers(Accept = 'application/json') |> 
    httr2::req_url_query(page_size = 1000) |> #size = 1000 so no need for pagination
    httr2::req_perform()
  
  # store result
  if(httr2::resp_is_error(req)){ # error handling
    stop(paste("could not find names for IDs. Error number: ", resp_status()))
  } else{
    resps <- httr2::resp_body_json(req)
    taxonomy_information <- resps$taxonomy_nodes 
    # extract counts, cant guarantee that each field will always have 7
    all_count_types <- unique(unlist(lapply(taxonomy_information, function(x) {
      if (is.null(x$taxonomy$counts)) {
        return(NULL)
      } else {
      sapply(x$taxonomy$counts, function(count) {
        count$type
        })
      }
    })))
    # Extract data from each taxonomy
    result_list <- lapply(taxonomy_information, function(x) {
      tax_data <- x$taxonomy
      
      # Skip if no taxonomy data
      if (is.null(tax_data)) {
        return(NULL)
      }
      # extract easy fields
      top_level <- data.frame(
        tax_id = tax_data$tax_id %||% NA,
        organism_name = tax_data$organism_name %||% NA,
        blast_name = tax_data$blast_name %||% NA,
        rank = tax_data$rank %||% NA,
        genomic_moltype = tax_data$genomic_moltype %||% NA,
        stringsAsFactors = FALSE
      )
      
      # create intermediary frame for count information
      count_list <- setNames(as.list(rep(0, length(all_count_types))), all_count_types)
      counts_df <- data.frame(count_list)
      
      # Fill in actual count values
      if (!is.null(tax_data$counts)) {
        for (count_entry in tax_data$counts) {
          counts_df[[count_entry$type]] <- count_entry$count
        }
      }
      
      # Combine everything
      result_df <- cbind(top_level, counts_df)
      
      # Add columns that are going to remain lists
      result_df$lineage <- list(tax_data$lineage)
      result_df$children <- list(tax_data$children)
      
      return(result_df)
    })
    
    # Remove NULL entries and bind
    result_list <- result_list[!sapply(result_list, is.null)]
    final_frame <- do.call(rbind, result_list)
    
    saveRDS(final_frame, file = paste0(file_location, "taxonomy_information_", 
      paste0(final_frame$organism_name, collapse = "_"), "_", Sys.Date(), ".rds"))
  }
  # file_location <- "../01_inputs/03_metadata/"
}
  


# 03. get metadata by taxon ID --------------------------------------------

# e.g. I run the first function to get the tax ID, then I can use that
# to lookup all accessions related to that ID (i.e. all Sphingomonas accessions)

# another NCBI REST API, with different params in the "url_path_append"

download_accession_metadata_by_taxID <- function(tax_ID, file_location){
  # throw an error if "taxonomic_name" is not a character
  stopifnot(
    "taxonomic_ID must be integer" = is.numeric(tax_ID),
    "taxonomic_ID cannot be empty" = length(tax_ID) > 0,
    "file_location cannot be empty" = length(file_location) > 0)
  
  # check directory "file_location" exists, if not, create
  source("../00_scripts/R_scripts/Utility_functions.R")
  check_and_create_dir(file_location)
  ncbi_secrets <- get_secrets("ncbi")
  
  # create hollow list
  results_list <- list()
 
  for (taxon in tax_ID) {
    # taxon <- 43668 / file_path <- "../01_inputs/03_metadata_copy/"
    # tax_ID <- c(43668, 1696)
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
     # path =  paste0(file_location, taxon, "_{Sys.Date()}_page{i}.json")
    )
    response_bodies <- purrr::map(resps, httr2::resp_body_json)
    # Extract all reports and combine into one list
    all_reports <- unlist(lapply(response_bodies, function(x) x$reports),
                          recursive = FALSE)
    names(all_reports) <- sapply(all_reports, function(x) x$accession)
    
    results_list[[as.character(taxon)]] <- all_reports
    message(paste(taxon, " has been processed"))
  }
  message("all done successfully")
  
  saveRDS(results_list, 
          file = paste0(file_location, "accession_information_", 
                        paste0(tax_ID, collapse = "_"),
                        "_", Sys.Date(), ".rds"))
}


# 04. get protein and/or genome sequence using accession ID ---------------------------------------------------


get_dataset_by_accession <- function(accession, ftype = "both"){
  # debug: accession = "GCF_002407065.1"
  # throw an error if "accession" is not a character
  stopifnot(
    "accession must be character" = is.character(accession),
    "accession cannot be empty" = length(accession) > 0,
    "type not in list of compatible types" = ftype %in% 
      c("genome", "protein", "both"))
  
  # NOTE: source path will have to be adjusted when moved to hawk
  source("../00_scripts/R_scripts/Utility_functions.R")
  ncbi_secrets <- get_secrets("ncbi")
  
  # convert "type" to formats expected by the API
  types <- data.frame(type = c("genome", "protein"), 
                      NCBI_value = c("GENOME_FASTA", "PROT_FASTA"),
                      extension = c("fna", "faa"))
  
  

  # setup list for bringing down fastas
  fasta_list <- list(protein = list(),
                     genome = list())
  
  if(ftype == "both"){
    fasta_type <- paste0(types$NCBI_value, collapse = ",")
    search_string <- paste0(types$extension, collapse = "|")
  } else{
    fasta_type <- types[types$type == ftype, "NCBI_value"]
    search_string <- types[types$type == ftype, "extension"]
    fasta_list <- fasta_list[ftype]
  } 
  
  
  # setup temp file for temporarily storing zip right after call
  temp_zip <- tempfaccessiontemp_zip <- tempfile(fileext = ".zip")
  
  # query NCBI database using REST API
  req <- httr2::request(ncbi_secrets$host) |> 
    httr2::req_auth_bearer_token(ncbi_secrets$api_key) |> 
    httr2::req_url_path_append('genome', 'accession', accession, 'download') |> 
    httr2::req_headers(Accept = 'application/zip') |> 
    httr2::req_url_query(include_annotation_type = fasta_type) |> 
    httr2::req_perform(path = temp_zip)
  
  # Find sequence files, read them into a frame
  all_files <- unzip(temp_zip, list = TRUE)$Name
  fasta_files <- data.frame(filename = all_files[grepl(paste0("\\.(", search_string, ")$"), all_files)])
  
  # error handling
  if(nrow(fasta_files) < 1){
    stop("No fasta files returned from API call")
  }
  
  # merge filename frame with "types"
  fasta_files$extension <- tools::file_ext(fasta_files$filename)
  types <- merge(types, fasta_files, by = "extension")
  
  
  # read in files matching "types" from "temp_zip" into "fasta_list"
  for(type in types$type){
   fasta_content <- readLines(unz(temp_zip, types[types$type == type, "filename"]))
   fasta_list[[type]] <- fasta_content
  }
  unlink(temp_zip)
  return(fasta_list)
}

