
# libraries ---------------------------------------------------------------


library(httr2)      # used to access NCBI REST API
library(tidyverse)  # all purpose tool
library(tidyjson)   # tidyverse extension for .json files
library(jsonlite)   # read in .json files
library(RPostgres)  # allows issueing of commands to postgres database
library(DBI)        # controls connection to database

# setup -------------------------------------------------------------------

# load in secrets (passwords etc) needed to link to stuff
secrets <- read.delim("../98_config/privatestuff.csv", header = T, sep = ",")
urlstring <- 'https://api.ncbi.nlm.nih.gov/datasets/v2/'
api_header <- secrets %>%
  filter(username == "ncbi") %>% pull(password) %>% as.character()
database_password <- secrets %>%
  filter(username == "postgres") %>% pull(password) %>% as.character()

# create frame of genus names
# the numbers are the NCBI IDs for the genera, which I looked up manually
# I now know it is possible to find them with the API
existing_taxa <- data.frame(taxon_id = c(53335, 13687, 1696, 43668, 33882), 
                            genus = c("pantoea", "sphingomonas", "brevibacterium", "brachybacterium",
                                      "microbacterium"))
# chocolatum was added for testing as there are many fewer samples, also funny name

# 01. Query the NCBI API --------------------------------------------------


for (i in 1:nrow(existing_taxa)) {
 # i <- 1 # debug line
  taxon <- existing_taxa$taxon_id[i] # set variables for URL
  genus <- existing_taxa$genus[i] # set variable to create file names
  
  req <- request(urlstring) %>% # link to NCBI database
    req_auth_bearer_token(api_header) %>% # prove I am trusted by site
    req_url_path_append('genome', 'taxon', taxon, 'dataset_report') %>% # make the command
    req_url_query(
      page_size = 200, # return 200 results per request, so I dont breach max query size
      filters.exclude_paired_reports = TRUE) %>% # bug fix, no longer bring down doubles (refer to week_4.qmd)
    req_headers(Accept = 'application/json') # take the output in .json format
  # now request is forged, but hasnt been sent, thus, the next block
  # taken from httr2 manual on how to do iterative requests
  
  resps <- req_perform_iterative( # will do multiple pages (size set in prev block)
    req,  
    next_req = iterate_with_cursor( 
      'page_token', 
      resp_param = function(resp) {
        resp_body_json(resp)$next_page_token
      } # method supported by httr2 to ask for next page of data, just copied from manual
    ),
    path = paste0("../02_analysis/01_existing_genera/", genus, "_", "{Sys.Date()}_page{i}.json")  # saves output to file
  )
}
# files will be moved to SQL, but if I want to reload later, then going straight
# would mean having to rerun the whole call again, so saving to files right thing to do

# 02. load output files to SQL database -----------------------------------

# create connection to database and create table if doesnt exist already
conn <- dbConnect(RPostgres::Postgres(), # begin connection
                  dbname = 'postgres', # name I made for database
                  host = 'localhost', # default host name (meaning my computer)
                  port = 5432, # default port for database
                  user = 'postgres', # root user 
                  password = database_password) # password

if (!dbExistsTable(conn, "genome_reports")) { # looks for specific table in the database
  dbExecute(conn, paste0( # if it doesnt exist then run this SQL command
    "CREATE TABLE ", "genome_reports", " (
        id SERIAL PRIMARY KEY,
        accession TEXT UNIQUE NOT NULL,
        genus_id TEXT NULL,
        report_data JSONB NOT NULL,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
      )"
  )) # This creates the table with all the collumns and the necessary charicteristics
} # this will mean the command wont throw an error

# Initialise a progress bar
pb <- txtProgressBar(min = 0, max = nrow(existing_taxa), style = 3)

# table now definitely exists, just have to populate it
for (i in 1:nrow(existing_taxa)) {
  # i <- 1 # debug line
  taxon <- existing_taxa$taxon_id[i] # set variables for URL
  genus <- existing_taxa$genus[i] # set variable to create file names
  
  file_pattern <- paste0(genus, "_*.json") # look for files starting with spec genus name, then _(anything).json ...
  directory <- "../02_analysis/01_existing_genera/" # ... in this directory
  # Get all file paths matching the pattern
  file_paths <- list.files(directory, pattern=glob2rx(file_pattern), full.names = TRUE) # glob2rx helps with REGEX, reads the names into a vector
  json_list <- lapply(file_paths, jsonlite::read_json) # takes names, reads in associated .json files to a list
  # json_list contains 1 entry per file and each file contains up to 200 samples. So, there is a list within a list
  # I want to write the samples individually into the database, the middle for loop loops over each file and pulls
  # out the reports object, which contains the sample datasets, the inner for loop reads each sample in turn and stores it in Db
  for (k in seq_along(json_list)) {
    json_reports <- json_list[[k]]$reports
    for (j in seq_along(json_reports))
    {
      
      single_report <- json_reports[[j]]
      report_json <- jsonlite::toJSON(single_report, auto_unbox = TRUE, pretty = TRUE) # reads indv things to a json objects, auto_unbox gives simpler str
      accession <- single_report$accession # accession pulled out because gets its own column
      
      dbExecute(
        conn,
        paste0("INSERT INTO ", "genome_reports", " (accession, genus_id, report_data) VALUES ($1, $2, $3) 
              ON CONFLICT (accession) DO UPDATE SET report_data = $3"), # if already exists, will update
        params = list(accession, taxon, report_json))
      
    }
  }
  # Print progress
  setTxtProgressBar(pb, i)
  
}

close(pb)

# 03. useful extra data --------------------------------------------------


existing_taxa_vector <- paste(existing_taxa$taxon_id, collapse = ",") # this call needs special formatting

req <- request(urlstring) %>% # link to NCBI database
    req_auth_bearer_token(api_header) %>% # prove I am trusted by site
    req_url_path_append('taxonomy', 'taxon', existing_taxa_vector, 'dataset_report')  %>% 
    req_headers(Accept = 'application/json') %>%  # take the output in .json format
  req_perform(path = paste0("../02_analysis/01_existing_genera/", Sys.Date(), "_existing_taxa.json"))
  # now request is forged, but hasnt been sent, thus, the next block
  # taken from httr2 manual on how to do iterative requests

if (!dbExistsTable(conn, "genus_info")) { # looks for specific table in the database
  dbExecute(conn, paste0( # if it doesnt exist then run this SQL command
    "CREATE TABLE ", "genus_info", " (
        id SERIAL PRIMARY KEY,
        genus_id TEXT UNIQUE NOT NULL,
        genus_name TEXT NULL,
        report_data JSONB NOT NULL,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
      )"
  )) # This creates the table with all the columns and the necessary characteristics
} # this will mean the command wont throw an error

json_list <- jsonlite::read_json(paste0("../02_analysis/01_existing_genera/", Sys.Date(), "_existing_taxa.json"))
pb <- txtProgressBar(min = 0, max = nrow(existing_taxa), style = 3)

json_reports <- json_list$reports
for (k in seq_along(json_reports)) {
  single_report <- json_reports[[k]]$taxonomy
  report_json <- jsonlite::toJSON(single_report, auto_unbox = TRUE, pretty = TRUE)
  genus_id <- single_report$tax_id
  genus_name <- single_report$current_scientific_name$name
    
  dbExecute(
    conn,
    paste0("INSERT INTO ", "genus_info", " (genus_id, genus_name, report_data) VALUES ($1, $2, $3) 
              ON CONFLICT (genus_id) DO UPDATE SET report_data = $3"), # if already exists, will update
    params = list(genus_id, genus_name, report_json))
    
  # Print progress
  setTxtProgressBar(pb, k)
}
close(pb)
