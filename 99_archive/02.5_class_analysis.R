# libraries ---------------------------------------------------------------


library(tidyverse)  # all purpose tool
library(RPostgres)  # allows issueing of commands to postgres database
library(DBI)        # controls connection to database
library(gt)         # more flexible tabling
library(rentrez)    # interacts with older NCBI API, connects hosts to classes
library(jsonlite)   # .json reading
library(XML)        # rentrez brings back .nxml, need this to read properly
library(httr2)      # NCBI REST API

# setup -------------------------------------------------------------------

# load in secrets (passwords etc) needed to link to stuff
secrets <- read.delim("../98_config/privatestuff.csv", header = T, sep = ",")
urlstring <- 'https://api.ncbi.nlm.nih.gov/datasets/v2/'
api_header <- secrets %>%
  filter(username == "ncbi") %>% pull(password) %>% as.character()
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


# 05. manual work on duff values ------------------------------------------


# I went through the 180 duff values manually looking for things to remove
dbExecute(conn,
"DELETE FROM public.unprocessed_hosts where host_value IN ('missing', 'not applicable', 'not determined', 'na', 'not available: not collected', 
                                                   'envirnment', 'environment', 'environmental', 'not provided', 'water from a lake', 'ucc strain', 'biological soil crusts', 'biofilm', 
                                                   'grown in blood agar culture medium', 'grown in blood agar culture medium.', 'host', 'natural / free-living', 
                                                   'sw4848was obtained from city park lake, baton rouge, louisiana, u.s.a.', 'zymobiomics microbial community standard strain',
                                                   'atcc strain', 'food', 'animal', 'laboratory', 'land crag', 'milk', 'nist mixed microbial rm strain', 'sediment of yellow river', 
                                                   'sewage', 'wildlife', 'digitaria sp.', 'freshwater lake', 'haemophisalus','living seaweed', 'not collected', 'plant', 'rain water',
                                                   'soil', 'sponge');")
# then I went through and set the alias for many human varients to "homo sapiens"

dbExecute(conn, 
"UPDATE public.unprocessed_hosts
SET alias = 'homo sapiens', host_id = '9606'
WHERE host_value LIKE '%homo%' OR host_value LIKE '%human%' OR host_value LIKE '%patient%'
OR host_value = 'soldier';")

# I did some testing on the NCBI REST API and found that more of my values could
# be done automatically. However, there are a few I would like to do manually first
# such as "bat", which is known to be difficult

dbExecute(conn, 
          "UPDATE public.unprocessed_hosts
SET alias = 'chiroptera', host_id = '9397'
WHERE host_value LIKE '%bat%';")

# now I need a dataframe for those without a host ID
result <- dbGetQuery(conn, 
                     "SELECT * FROM public.unprocessed_hosts WHERE host_id IS NULL;")

# manually going down list
host_list <- arrange(result, host_value) %>% 
  select(host_value)
host_list

# now I need to store the results
x <- data.frame(paste0("'", host_list$host_value, "', "))
print(x, row.names = FALSE)

# list was created manually in text file "02_analysis/03_unknown_hosts/manually_assigned_hosts.txt"

manual_hosts <- read.table("../02_analysis/03_unknown_hosts/manually_assigned_hosts.txt", 
                 sep = ",", 
                 quote = "'", 
                 col.names = c("host", "code"),
                 stringsAsFactors = FALSE)
manual_hosts$host <- trimws(manual_hosts$host)
manual_hosts$code <- trimws(manual_hosts$code)

dbWriteTable(conn, "temp_hosts_codes", manual_hosts, temporary = TRUE, overwrite = TRUE)
# bring two tables together in database
res <- dbGetQuery(conn, "SELECT * FROM temp_hosts_codes tc
                         INNER JOIN unprocessed_hosts uh
                         ON uh.host_value = tc.host;")

dbExecute(conn, "UPDATE unprocessed_hosts
                 SET host_id = code
                 FROM temp_hosts_codes 
                 WHERE host_value = host;")


# 06. looking up the taxonomy information ---------------------------------


res <- dbGetQuery(conn, "SELECT * FROM  unprocessed_hosts;")

pb <- txtProgressBar(min = 0, max = nrow(res), style = 3)

for (i in 1:nrow(res)) {
   #i <- 1 
# retrieve tax record (also from rentrez manual)
  taxid <- res$host_id[i]
  host_value <- res$host_value[i]
 
  tax_rec <- entrez_fetch(db="taxonomy", id=taxid, rettype="xml", parsed=TRUE)
    
  tax_list <- XML::xmlToList(tax_rec)
  tax_json <- jsonlite::toJSON(tax_list, auto_unbox = TRUE)
  
  dbExecute(
    conn,
    paste0("INSERT INTO ", "host_reports", " (host_value, host_id, report_data) VALUES ($1, $2, $3) 
              ON CONFLICT (host_value) DO UPDATE SET report_data = $3"), # if already exists, will update
    params = list(host_value, taxid, tax_json))
  
  # Print progress
  setTxtProgressBar(pb, i)
  Sys.sleep(0.05)
}

close(pb)


# 07. classification ------------------------------------------------------

# now all the duff samples have been turned good and added to the list, I am ready
# to create my groups.

# gets every host and its taxonomic class
host_to_class_info <- dbGetQuery(conn, 
      "SELECT host_value,
        report_data ->'Taxon'->>'ScientificName' as taxon_name,
        report_data ->'Taxon'->>'Rank' as taxon_rank,
        (value->>'ScientificName') as class_name,
        (value->>'TaxId') as class_tax_id
        FROM host_reports,
        LATERAL jsonb_each(report_data->'Taxon'->'LineageEx') as t(key, value)
        WHERE value->>'Rank' = 'class';")

# gets every accession with a valid host and the host name
accession_to_host_info <- dbGetQuery(conn, "SELECT Q1.genus_id, Q1.genus_name, Q1.accession, Q1.host FROM 
(SELECT gi.genus_id, gi.genus_name, gr.accession, gr.report_data #>> '{assembly_info, biosample, host}' as host
  FROM public.genome_reports gr
  INNER JOIN public.genus_info gi ON gi.genus_id = gr.genus_id) AS Q1
WHERE Q1.host IS NOT NULL AND
Q1.host NOT IN ('missing', 'mssing', 'missed', 'not applicable', 'unknown', 'N/A', 'not collected', 'no collected', 'not available', 'Not applicable');")

# make all words lower case so that I dont get "Bat" and "bat" as seperate entries
accession_to_host_info <- accession_to_host_info %>% 
  mutate(host_value = tolower(host))

# so this is my grand-total object of all the important data
accession_to_class_info <- host_to_class_info %>% 
  inner_join(accession_to_host_info, by= join_by(host_value))

# now I can see how many accessions there are of each Class
accession_by_class_count <- accession_to_class_info %>%  
  count(class_name, sort = T)

# Class count, but by bacterial genus
accession_by_class_and_genus_count <- accession_to_class_info %>% 
  group_by(class_name, genus_name) %>% 
  count() %>% 
  arrange(class_name, desc(n))

# now pivoted
pivoted_accession_by_class_and_genus_count <- accession_by_class_and_genus_count %>% 
  pivot_wider(names_from = genus_name, values_from = n)

# Class count, but by bacterial genus, prob best for display purposes
accession_by_class_and_genus_count_for_display <- accession_to_class_info %>% 
  group_by(class_name, genus_name) %>% 
  count(name = "count_of_accessions") %>% 
  arrange(desc(count_of_accessions))

# now this pivoted
pivoted_accession_by_class_and_genus_count_for_display <- accession_by_class_and_genus_count_for_display %>% 
  pivot_wider(names_from = genus_name, values_from = count_of_accessions)

# TODO: write out accession_to_class_info to txt file so I can do these tables in notebook
write_delim(accession_to_class_info, "../03_outputs/host_evaluation/accession_to_host_class.tsv", delim = "\t")

# 99. legacy code for using REST API --------------------------------------

# wasnt very reliable in testing, just going to do manually

# now I will use the REST API to see how many more I can automatically fill in.
# # First, I need to pull the list from SQL
# result <- dbGetQuery(conn, 
#                      "SELECT * FROM public.unprocessed_hosts WHERE host_id IS NULL;")
# 
# dbExecute(conn,
#           "ALTER TABLE public.unprocessed_hosts ADD returned_rank VARCHAR(100);")
# 
# pb <- txtProgressBar(min = 0, max = nrow(result), style = 3)
# # then run it
# for (i in 1:nrow(result)) {
#   i <- 2
#   taxon_query <- result$host_value[i]
#   taxon_query_url <- gsub(" ", "%20", taxon_query)
#   req <- request(urlstring) %>% # link to NCBI database
#   req_auth_bearer_token(api_header) %>% # prove I am trusted by site
#   req_url_path_append('taxonomy', 'taxon_suggest', taxon_query_url)  %>%
#   req_headers(Accept = 'application/json') %>%  # take the output in .json format
#   req_perform()
#   resp <- resp_body_json(req)
#   
#   if (length(resp) > 0) {
#     filtered_resp <- map_dfr(resp$sci_name_and_ids, ~ as_tibble(.x)) %>%
#       filter(!is.na(rank), rank != "") %>%
#       slice(1)
#     
#     host_name <- filtered_resp$sci_name
#     host_id <- filtered_resp$tax_id
#     returned_rank <- filtered_resp$rank
#     
#     # update table
#     dbExecute(
#       conn,"UPDATE public.unprocessed_hosts 
#             SET alias = $1, host_id = $2, returned_rank = $3
#             WHERE host_value = $4;",
#       params = list(host_name, host_id, returned_rank, taxon_query))
#   }
#   # Print progress
#   setTxtProgressBar(pb, i)
# }
# 
# close(pb)
