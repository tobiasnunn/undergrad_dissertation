# libraries ---------------------------------------------------------------


library(tidyverse)  # all purpose tool
library(RPostgres)  # allows issueing of commands to postgres database
library(DBI)        # controls connection to database
library(gt)         # more flexible tabling
library(rentrez)    # interacts with older NCBI API, connects hosts to classes
library(jsonlite)   # .json reading
library(XML)        # rentrez brings back .nxml, need this to read properly

# setup -------------------------------------------------------------------

# load in secrets (passwords etc) needed to link to stuff
secrets <- read.delim("../98_config/privatestuff.csv", header = T, sep = ",")

database_password <- secrets %>%
  filter(username == "postgres") %>% pull(password) %>% as.character()

# 01. bring in data ------------------------------------------------------

conn <- dbConnect(RPostgres::Postgres(), 
                  dbname = 'postgres',
                  host = 'localhost', 
                  port = 5432, 
                  user = 'postgres', 
                  password = database_password) 

# get data from database
result <- dbGetQuery(conn, "SELECT Q1.genus_id, Q1.genus_name, Q1.accession, Q1.host FROM 
(SELECT gi.genus_id, gi.genus_name, gr.accession, gr.report_data #>> '{assembly_info, biosample, host}' as host
  FROM public.genome_reports gr
  INNER JOIN public.genus_info gi ON gi.genus_id = gr.genus_id) AS Q1
WHERE Q1.host IS NOT NULL AND
Q1.host NOT IN ('missing', 'mssing', 'missed', 'not applicable', 'unknown', 'N/A', 'not collected', 'no collected', 'not available', 'Not applicable');")

# make all words lower case so that I dont get "Bat" and "bat" as seperate entries
lower <- result %>% 
  mutate(host_lower = tolower(host))

# unique the list
unique_host <- data.frame(host = unique(lower$host_lower), 
                          taxonomy_id = as.integer(NA),
                          division = NA,
                          taxonomy_data = NA) 


# 02. test querying ------------------------------------------------------

pb <- txtProgressBar(min = 0, max = nrow(unique_host), style = 3)

for (i in 1:nrow(unique_host)) {
  # i <- 3 debug line
  # test whether one word entries might work, specifically "eucalyptus" and "medicago"
  species_name <- unique_host[i,1] # name of species to pass onto API
  
  # search for taxonomy ID on NCBI website (from rentrez manual)
  search_result <- entrez_search(db = "taxonomy", term = species_name, retmax = 1)
  
  # retrieve tax record (also from rentrez manual)
  taxid <- search_result$ids[1]
  if (!is.null(taxid[[1]])) {
    tax_rec <- entrez_fetch(db="taxonomy", id=taxid, rettype="xml", parsed=TRUE)
    
    tax_list <- XML::xmlToList(tax_rec)
    tax_json <- jsonlite::toJSON(tax_list, auto_unbox = TRUE)
    tax_division <- tax_list$Taxon$Division
    
    unique_host$division[i] <- tax_division
    unique_host$taxonomy_data[i] <- tax_json
    unique_host$taxonomy_id[i] <- taxid
  }
  # Print progress
  setTxtProgressBar(pb, i)
  Sys.sleep(0.05)
}

close(pb)

# exploration
divisions <- unique(unique_host$division)
divisions
na_count <- sum(is.na(unique_host$division))
na_count
na_prop <- na_count / nrow(unique_host)
na_prop


# 03. write to SQL database -----------------------------------------------

# this table will just contain the ones that returned a value and will
# leave the ones that are NULL for a later table after manual work
if (!dbExistsTable(conn, "host_reports")) { 
  dbExecute(conn, paste0( 
    "CREATE TABLE ", "host_reports", " (
        id SERIAL PRIMARY KEY,
        host_value TEXT UNIQUE NOT NULL,
        host_id TEXT NOT NULL,
        report_data JSONB NOT NULL,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
      )"
  )) 
} 

# now populate table
filtered_unique_host <- filter(unique_host, !is.na(taxonomy_id))
pb <- txtProgressBar(min = 0, max = nrow(filtered_unique_host), style = 3)

for (i in 1:nrow(filtered_unique_host)) {
  #i <- 1
  single_report <- filtered_unique_host$taxonomy_data[i]
  host <- filtered_unique_host$host[i]
  taxon_id <- filtered_unique_host$taxonomy_id[i]
  
  dbExecute(
    conn,
    paste0("INSERT INTO ", "host_reports", " (host_value, host_id, report_data) VALUES ($1, $2, $3) 
              ON CONFLICT (host_value) DO UPDATE SET report_data = $3"), # if already exists, will update
    params = list(host, taxon_id, single_report))
  # Print progress
  setTxtProgressBar(pb, i)
}

close(pb)


# 04. table of NULL host_reports ------------------------------------------


# this table will just contain the ones that returned a value and will
# leave the ones that are NULL for a later table after manual work
if (!dbExistsTable(conn, "unprocessed_hosts")) { 
  dbExecute(conn, paste0( 
    "CREATE TABLE ", "unprocessed_hosts", " (
        id SERIAL PRIMARY KEY,
        host_value TEXT UNIQUE NOT NULL,
        host_id TEXT NULL,
        alias TEXT NULL,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
      )"
  )) 
} 

# now populate table
filtered_unique_host <- filter(unique_host, is.na(taxonomy_id))
pb <- txtProgressBar(min = 0, max = nrow(filtered_unique_host), style = 3)

for (i in 1:nrow(filtered_unique_host)) {
  #i <- 1
  host <- filtered_unique_host$host[i]
  
  dbExecute(
    conn,
    paste0("INSERT INTO ", "unprocessed_hosts", " (host_value) VALUES ($1)"), 
    params = list(host))
  # Print progress
  setTxtProgressBar(pb, i)
}

close(pb)

