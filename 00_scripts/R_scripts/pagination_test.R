download_taxonomic_information_by_tax_ID <- function(tax_ID, file_location){
  # [Your existing validation code stays the same]
  stopifnot(
    "taxonomic_ID must be integer" = is.numeric(tax_ID),
    "taxonomic_ID cannot be empty" = length(tax_ID) > 0,
    "file_location cannot be empty" = length(file_location) > 0)
  
  # [Your existing setup stays the same]
  source(here::here("00_scripts/R_scripts/Utility_functions.R"))
  check_and_create_dir(file_location)
  ncbi_secrets <- get_secrets("ncbi")
  tax_ID_string <- paste0(tax_ID, collapse = ",")
  
  # NEW: Initialize pagination variables
  all_taxonomy_information <- list()  # Store all results
  next_page_token <- NULL              # Start with no token
  page_number <- 1                     # Track which page we're on
  
  # NEW: Loop until no more pages
  repeat {
    cat("Fetching page", page_number, "...\n")
    
    # Build request (mostly your existing code)
    req <- httr2::request(ncbi_secrets$host) |> 
      httr2::req_auth_bearer_token(ncbi_secrets$api_key) |> 
      httr2::req_url_path_append('taxonomy', 'taxon', tax_ID_string) |> 
      httr2::req_headers(Accept = 'application/json') |> 
      httr2::req_url_query(page_size = 1000)
    
    # NEW: Add page_token if we have one (skip on first iteration)
    if (!is.null(next_page_token)) {
      req <- req |> httr2::req_url_query(page_token = next_page_token)
    }
    
    # Perform request
    req <- httr2::req_perform(req)
    
    # [Your existing error handling]
    if(httr2::resp_is_error(req)){
      stop(paste("Could not find names for IDs. Error number:", httr2::resp_status(req)))
    }
    
    # Parse response
    resps <- httr2::resp_body_json(req)
    
    # NEW: Check if we got any results
    if (is.null(resps$taxonomy_nodes) || length(resps$taxonomy_nodes) == 0) {
      cat("No more results found\n")
      break
    }
    
    # Store this page's results
    all_taxonomy_information[[page_number]] <- resps$taxonomy_nodes
    cat("Retrieved", length(resps$taxonomy_nodes), "records\n")
    
    # NEW: Check for next page token
    next_page_token <- resps$next_page_token
    
    if (is.null(next_page_token) || next_page_token == "") {
      cat("No more pages available\n")
      break
    }
    
    page_number <- page_number + 1
    
    # Be nice to the API
    Sys.sleep(0.35)
  }
  
  # NEW: Combine all pages
  taxonomy_information <- unlist(all_taxonomy_information, recursive = FALSE)
  
  cat("Total records retrieved:", length(taxonomy_information), "\n")
  
  # [Rest of your existing processing code stays EXACTLY the same]
  all_count_types <- unique(unlist(lapply(taxonomy_information, function(x) {
    if (is.null(x$taxonomy$counts)) {
      return(NULL)
    } else {
      sapply(x$taxonomy$counts, function(count) {
        count$type
      })
    }
  })))
  
  result_list <- lapply(taxonomy_information, function(x) {
    tax_data <- x$taxonomy
    
    if (is.null(tax_data)) {
      return(NULL)
    }
    
    top_level <- data.frame(
      tax_id = tax_data$tax_id %||% NA,
      organism_name = tax_data$organism_name %||% NA,
      blast_name = tax_data$blast_name %||% NA,
      rank = tax_data$rank %||% NA,
      genomic_moltype = tax_data$genomic_moltype %||% NA,
      stringsAsFactors = FALSE
    )
    
    count_list <- setNames(as.list(rep(0, length(all_count_types))), all_count_types)
    counts_df <- data.frame(count_list)
    
    if (!is.null(tax_data$counts)) {
      for (count_entry in tax_data$counts) {
        counts_df[[count_entry$type]] <- count_entry$count
      }
    }
    
    result_df <- cbind(top_level, counts_df)
    result_df$lineage <- list(tax_data$lineage)
    result_df$children <- list(tax_data$children)
    
    return(result_df)
  })
  
  result_list <- result_list[!sapply(result_list, is.null)]
  final_frame <- do.call(rbind, result_list)
  
  saveRDS(final_frame, file = paste0(file_location, "taxonomy_information_", 
                                     paste0(final_frame$organism_name, collapse = "_"), 
                                     "_", Sys.Date(), ".rds"))
  
  return(final_frame)  # NEW: Return the data so you can check it
}