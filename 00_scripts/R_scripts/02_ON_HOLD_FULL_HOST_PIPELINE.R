
# Introduction ------------------------------------------------------------


# the preceding 7 files in the under the banner "02" were used for development
# I plan to use this file to collate all the steps into a more clear pipeline
# for when I will want to run the steps in the future


# libraries ---------------------------------------------------------------

library(tidyverse)  # all purpose tool
library(RPostgres)  # allows issueing of commands to postgres database
library(DBI)        # controls connection to database

# setup -------------------------------------------------------------------

# I bring in a lot of data from a file called /98_config/privatestuff_v2.csv
# this makes the setup cleaner and protects some sensative data of mine
# the file will not appear in the repo for this reason, basically
# I have a column for every field in the setup section here.

secrets <- read.delim("../98_config/privatestuff_v2.csv", header = T, sep = ",")

#separate secrets for different sections
postgres_secrets <- filter(secrets, purpose == "database")
ncbi_secrets <- filter(secrets, purpose == "ncbi")

# forge connection with database
conn <- dbConnect(RPostgres::Postgres(), # begin connection
                  dbname = postgres_secrets$dbname, # name I made for database
                  host = postgres_secrets$host, # default host name (meaning my computer)
                  port = postgres_secrets$port, # default port for database
                  user = postgres_secrets$username, # root user 
                  password = postgres_secrets$password, # password
                  base::list(sslmode="prefer", connect_timeout="10"), 
                  service = NULL) # secures connection to database, good for remote databases

# cleanup
rm(secrets, postgres_secrets)
# im leaving ncbi_secrets as I will need it multiple times later


# 01. initial provision of genus information ------------------------------

# create frame of genus names
# the numbers are the NCBI IDs for the genera, which I looked up manually
# I now know it is possible to find them with the API
taxon_id <- c(53335, 13687, 1696, 43668, 33882) #got from NCBI Taxonomy browser
taxon_name <- c("Pantoea", "Sphingomonas", "Brevibacterium", "Brachybacterium",
                "Microbacterium") # as said by GTDB-Tk de_novo_wf last year
taxon_type <- c("genus") # at some point I may want to be working on families (maybe even species)
metadata_retrieved <- TRUE  # if running new groups, set to FALSE
bac_tax_data <- data.frame(taxon_id = taxon_id, 
                         taxon_name = taxon_name,
                         taxon_type = taxon_type,
                         metadata_retrieved = metadata_retrieved)
#clean
rm(taxon_id, taxon_name, taxon_type, metadata_retrieved)


# now I want to put that information into a table in the SQL database

#table creation
if (!dbExistsTable(conn, "taxa_of_study")) { # looks for specific table in the database
  dbExecute(conn, paste0( # if it doesnt exist then run this SQL command
    "CREATE TABLE ", "taxa_of_study", " (
        id SERIAL PRIMARY KEY,
        taxon_id INTEGER UNIQUE NOT NULL,
        taxon_name TEXT NOT NULL,
        taxon_type TEXT NOT NULL,
        metadata_retrieved BOOLEAN NOT NULL,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
      )"
  )) # This creates the table with all the collumns and the necessary charicteristics
} # this will mean the command wont throw an error

dbExecute(conn, "DELETE FROM taxa_of_study") # makes sure the table is empty before processing
# table population
dbWriteTable(conn, "taxa_of_study", bac_tax_data, append = T)
rm(bac_tax_data)

# 02. get metadata from NCBI API ------------------------------------------

genera_not_held <- dbGetQuery(conn, "SELECT * FROM public.taxa_of_study
                                     WHERE metadata_retrieved = FALSE;") 
for (i in 1:nrow(genera_not_held)) {
  # i <- 1 # debug line
  taxon_id <- genera_not_held$taxon_id[i] # set variables for URL
  taxon_name <- genera_not_held$taxon_name[i] # set variable to create file names
  
  req <- request(ncbi_secrets$host) %>% # link to NCBI database
    req_auth_bearer_token(ncbi_secrets$api_key) %>% # prove I am trusted by site
    req_url_path_append('genome', 'taxon', taxon_id, 'dataset_report') %>% # make the command
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

