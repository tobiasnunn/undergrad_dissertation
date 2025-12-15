# for both pseudomonas aer and pseudomonas putrida, get meta for the accessions that meet the following criteria:
# - complete genome
# - exclude atypical
# - 2015 or newer

# then pull out quality metrics - N50, number of contigs, contamination level so that we can select a reasonable number
# from each species

# get some basic metadata about the species
# get the reference samples for the spacies so we have the "total sequence length" to compare with
# get recent complete genomes for the species
# extract key quality metrics (completeness,contamination, N50, comp to reference length)
# order list by most desirable for each of those elements
# take the top 10-50 from each of those groups

library(httr2)
library(jsonlite)
library(tidyverse)

source(here::here("00_scripts/R_scripts/Utility_functions.R"))
ncbi_secrets <- get_secrets("ncbi")

species <- c("287", "303")
species_list <- paste(species, collapse = ',')
# get the taxon info for later


# get genus-level metadata ------------------------------------------------

req <- httr2::request(ncbi_secrets$host) |> 
  httr2::req_auth_bearer_token(ncbi_secrets$api_key) |>
  httr2::req_url_path_append('taxonomy', 'taxon', species_list,'dataset_report') |>
  httr2::req_headers(Accept = 'application/json') |>
  httr2::req_perform()



species_info_df <- map(resp_body_json(req)$reports, function(report) {
  tax <- report$taxonomy
  
  tibble(
    tax_id = tax$tax_id,
    rank = tax$rank,
    name = tax$current_scientific_name$name
  )
}) |> list_rbind()



# fetch the reference results ---------------------------------------------


results_list_reference <- list()

for(taxon in species){
  req <- httr2::request(ncbi_secrets$host) |> 
    httr2::req_auth_bearer_token(ncbi_secrets$api_key) |>
    httr2::req_url_path_append('genome', 'taxon', taxon, 'dataset_report') |>
    httr2::req_url_query(
      page_size = 200,
      filters.exclude_paired_reports = TRUE,
      filters.exclude_atypical = TRUE,
      filters.assembly_version = "current",
      filters.reference_only = TRUE) |>
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
  
  results_list_reference[[as.character(taxon)]] <- all_reports
}

# turn into dataframe
reference_df <- imap(results_list_reference, function(taxon_data, tax_id) {
  imap(taxon_data, function(accession_data, accession_id) {
    tibble(
      tax_id = as.integer(tax_id),
      accession_id = accession_id,
      total_sequence_length = pluck(accession_data, "assembly_stats", "total_sequence_length", .default = NA),
    )
  }) |> list_rbind()
}) |> list_rbind() |> select(-accession_id) |> 
  rename(total_reference_sequence_length = total_sequence_length)


# get recent complete genomes for 2 species -------------------------------


# 287 = aeruginosa, 303 = putida
results_list <- list()

for(taxon in species){
  req <- httr2::request(ncbi_secrets$host) |> 
    httr2::req_auth_bearer_token(ncbi_secrets$api_key) |>
    httr2::req_url_path_append('genome', 'taxon', taxon, 'dataset_report') |>
    httr2::req_url_query(
      page_size = 200,
      filters.exclude_paired_reports = TRUE,
      filters.exclude_atypical = TRUE,
      filters.assembly_version = "current",
      filters.assembly_level = "complete_genome",
      filters.first_release_date = "2015-01-01T00:00:00.000") |>
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
}


results_df <- imap(results_list, function(taxon_data, tax_id) {
  imap(taxon_data, function(accession_data, accession_id) {
    tibble(
      tax_id = as.integer(tax_id),
      accession_id = accession_id,
      number_of_contigs = pluck(accession_data, "assembly_stats", "number_of_contigs", .default = NA),
      contig_n50 = pluck(accession_data, "assembly_stats", "contig_n50", .default = NA),
      total_sequence_length = pluck(accession_data, "assembly_stats", "total_sequence_length", .default = NA),
      completeness = pluck(accession_data, "checkm_info", "completeness", .default = NA),
      contamination = pluck(accession_data, "checkm_info", "contamination", .default = NA),
    )
  }) |> list_rbind()
}) |> list_rbind() |> 
  left_join(species_info_df, by = join_by(tax_id)) |> 
  left_join(reference_df, by = join_by(tax_id)) |>
  select(name, everything()) |>
  mutate(length_comparison = as.numeric(total_sequence_length) / as.numeric(total_reference_sequence_length)) |>
  filter(completeness >= 95 & contamination <= 5 & number_of_contigs < 4) |>
  arrange(number_of_contigs, desc(contig_n50), desc(length_comparison))

summary_stat <- results_df |> count(name)
str(results_df)

# order list by most desirable for each of those elements
# take the top 10-50 from each of those groups

top50 <- results_df %>%
  group_by(name) %>%
  slice_head(n = 50) %>%
  ungroup()


# old stuff ---------------------------------------------------------------



# bind datarames together so get a full list
banana23 <- rbind(banana2, banana20)

bazinga <- banana23 |> count(name)
#write_csv(banana2, "banana.csv")
# Claude tells me that genomics studies take a representative 10-50 accessions from each group
# So, I will take the top 50 from each species in the banana2 object and bring in those accessions


# get reference genomes for comparing genome size
results_list_reference <- list()

for(taxon in species){
  req <- httr2::request('https://api.ncbi.nlm.nih.gov/datasets/v2/') |>
    httr2::req_auth_bearer_token(api_key) |>
    httr2::req_url_path_append('genome', 'taxon', taxon, 'dataset_report') |>
    httr2::req_url_query(
      page_size = 200,
      filters.exclude_paired_reports = TRUE,
      filters.exclude_atypical = TRUE,
      filters.assembly_version = "current",
      #filters.assembly_level = "contig",
      #filters.assembly_level = "complete_genome",
      #filters.first_release_date = "2015-01-01T00:00:00.000", 
      filters.reference_only = TRUE) |>
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
  
  results_list_reference[[as.character(taxon)]] <- all_reports
}

# completeness percentile
banana6 <- banana23 %>%
  group_by(name) %>%
  arrange(desc(completeness_percentile), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()
write_csv(banana6, "banana6.csv")