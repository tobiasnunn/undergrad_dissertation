
# 01. libraries -----------------------------------------------------------


library(httr2)
library(tidyverse)


# 02. obtain taxon IDs for the 5 genera -----------------------------------

# I will do this the fancy way. By using the NCBI REST API "taxon_suggest" endpoint
# e.g. by searching for "brevibacterium" and using a few parameters, 
# I can get the taxon ID to pass to the next call
# this will need to happen 5 times, as the API can only send one query at a time

# read in API key and URL using the function I just created
source("../00_scripts/R_scripts/Utility_functions.R") #bring in the functions script
# path relative to the .proj file

# call the function
ncbi_secrets <- get_secrets("ncbi")

# run the API call 5 times, I decided a for loop was best for this
# it will read all 5 outputs into 1 object

# vector of genus names
genera <- c("Brachybacterium", "Brevibacterium", "Microbacterium", "Pantoea", "Sphingomonas")

# create hollow for the loop
genus_info <- data.frame(tax_id = numeric(), sci_name = character())

# for loop
for (i in 1:length(genera)) {
  # connect to NCBI database, bring down data into "req" object
  req <- request(ncbi_secrets$host) %>% 
    req_auth_bearer_token(ncbi_secrets$api_key) %>% 
    req_url_path_append('taxonomy', 'taxon_suggest', genera[i]) %>% 
    req_headers(Accept = 'application/json') %>% 
    req_url_query(tax_rank_filter = "higher_taxon") %>% 
    req_perform()
  
  # store result
  result <- resp_body_json(req)$sci_name_and_ids %>%
    keep(~ .x$rank == "GENUS") %>% # using the purrr package to filter the list
    map_dfr(~ data.frame(tax_id = as.numeric(.x$tax_id), sci_name = .x$sci_name))
  
  genus_info <- rbind(genus_info, result) # bind the new row into the set_up dataframe
}


# 03. download metadata ---------------------------------------------------

# another NCBI REST API, with different params in the "url_path_append"
for (i in 1:nrow(genus_info)) {
  req <- request(ncbi_secrets$host) %>% 
    req_auth_bearer_token(ncbi_secrets$api_key) %>% 
    req_url_path_append('genome', 'taxon', genus_info$tax_id[i], 'dataset_report') %>% 
    req_url_query(
      page_size = 200, 
      filters.exclude_paired_reports = TRUE) %>% 
    req_headers(Accept = 'application/json')
  
  # keeps calling req until all pages are returned ("page_size  200")
  resps <- req_perform_iterative(
    req,  
    next_req = iterate_with_cursor(
      'page_token', 
      resp_param = function(resp) {
        resp_body_json(resp)$next_page_token
      }
    ),
    path = paste0("../01_inputs/03_metadata/", genus_info$sci_name[i], "_{Sys.Date()}_page{i}.json")  
  )
}
