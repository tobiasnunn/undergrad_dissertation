# try to fetch all 44000 pseudomonas aeriglaladlkj samples


# libraries ---------------------------------------------------------------


library(httr2)
library(jsonlite)
library(tidyverse)

api_key <- 'd4f921c91015c68f48f656a99ad5bcf84a08'
species <- c("287")
species_list <- paste(species, collapse = ',')


# taxon info --------------------------------------------------------------


# get the taxon info for later

req <- httr2::request('https://api.ncbi.nlm.nih.gov/datasets/v2/') |>
  httr2::req_auth_bearer_token(api_key) |>
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


# get 44k accessions for P. aeruginosa ------------------------------------

# done in batches of a few years to 1 year intervals so as to not hit the pull limit

# 287 = aeruginosa, 303 = putida
results_list <- list()

  req <- httr2::request('https://api.ncbi.nlm.nih.gov/datasets/v2/') |>
    httr2::req_auth_bearer_token(api_key) |>
    httr2::req_url_path_append('genome', 'taxon', taxon, 'dataset_report') |>
    httr2::req_url_query(
      page_size = 1000,
      filters.exclude_paired_reports = TRUE,
      filters.exclude_atypical = TRUE,
      filters.first_release_date = '2025-01-01T00:00:00.000',
      filters.last_release_date = '2025-12-31T00:00:00.000',
      filters.assembly_version = 'current') |>
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
  
  results_list[["2025"]] <- all_reports
  
  saveRDS(results_list, file = "results_list.rds")
  
  results_list[[as.character(taxon)]] <- all_reports


# find the host information -----------------------------------------------

# read in output from last section so no have to run again
# results_list <- readRDS("C:/Users/tobyn/source/undergrad_dissertation/01_inputs/03_metadata/results_list.rds")

banana2 <- imap(results_list, function(taxon_data, tax_id) {
  imap(taxon_data, function(accession_data, accession_id) {
    checkm_info <- pluck(accession_data, "checkm_info", .default = list())
    assembly_stats <- pluck(accession_data, "assembly_stats", .default = list())
    tibble(
      tax_id = pluck(accession_data, "organism", "tax_id", .default = NA),
      accession_id = accession_id,
      #number_of_contigs = pluck(accession_data, "assembly_stats", "number_of_contigs", .default = NA),
      host = pluck(accession_data, "assembly_info", "biosample", "host", .default = NA),
      iso_source = pluck(accession_data, "assembly_info", "biosample", "isolation_source", .default = NA)
    )
  }) |> list_rbind()
}) |> list_rbind()


count_of_thing <- banana2 %>% count(host)


# filter the list to just animal-hosts

# filter list to just animal hosts ----------------------------------------

# unique list to compare with
unique_list <- data.frame(unique(banana2$host))

#filtered_metadata <- filtered_metadata[-which(is.na(filtered_metadata$host)),]


#  method 1 - big list ----------------------------------------------------


# proper list to filter
filtered_metadata <- 
  banana2[which(banana2$host %in% 
                  c( #cow group
                    "Bos taurus", "Bos indicus", "Bos taurus indicus","bovine",
                    "cattle", "cow",
                    # dog group
                    "canis", "Canis", "canine", "Canine", "Canis lupus", 
                    "Canis lupus familiaris", "coyote",
                    # fish group
                    "Astronotus ocellatus", 
                    # goat group ( and sheep?)
                    "Capra hircus",
                    # cat group
                    "cat", "Cat",
                    # deer group
                    "Cervus elaphus",
                    # dolphin group
                    "Delphinidae",
                    # chicken group
                    "chicken",
                    # other-bird group (maybe duck if there are enough)
                    "Dendrocygna viduata")),]

# method 2 - many sub-frames that later get combined ----------------------


# cows
cow_group <- banana2[which(banana2$host %in% 
                  c("Bos taurus", "Bos indicus", "Bos taurus indicus","bovine",
                    "cattle", "cow")),] |> mutate(group = "Bovine_group")

# dogs
dog_group <- banana2[which(banana2$host %in% 
                             c("canis", "Canis", "canine", "Canine", "Canis lupus", 
                               "Canis lupus familiaris", "coyote")),] |>
  mutate(group = "Canine_group")

#cats
cat_group <- banana2[which(banana2$host %in% 
                             c("cat", "Cat")),] |>
  mutate(group = "Feline_group")

#chickens
chicken_group <- banana2[which(banana2$host %in% 
                                 c("chicken")),] |>
  mutate(group = "Chicken_group")

#fish
fish_group <- banana2[which(banana2$host %in% 
                                 c("Astronotus ocellatus")),] |>
  mutate(group = "Fish_group")

#goat
goat_group <- banana2[which(banana2$host %in% 
                              c("Capra hircus")),] |>
  mutate(group = "Goat_group")

#deer
deer_group <- banana2[which(banana2$host %in% 
                              c("Cervus elaphus")),] |>
  mutate(group = "Deer_group")

#dolphin
dolph_group <- banana2[which(banana2$host %in% 
                              c("Delphinidae")),] |>
  mutate(group = "Dolphin_group")

#duck
duck_group <- banana2[which(banana2$host %in% 
                               c("Dendrocygna viduata")),] |>
  mutate(group = "Duck_group")


# bind them all together

filtered_metadata2 <- rbind(cow_group, dog_group, cat_group, chicken_group,
                            fish_group, goat_group, deer_group, dolph_group,
                            duck_group)

# count by group
count_of_thing2 <- filtered_metadata2 %>% count(group)
