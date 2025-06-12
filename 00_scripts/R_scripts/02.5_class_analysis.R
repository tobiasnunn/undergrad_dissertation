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

'alhagi sparsifolia shap.', '218100'
'allium sativum var. ophioscorodon', '4682'
'allomyria dichotoma', '273928'
'alyssum serpyllifolium ssp. lusitanicum', '226043' 
'amphibian : african clawed toad', '8355'
'amphisylla', '184247'
'arabodopsis thaliana', '3702'
'badger', '9655'
'bannana', '4640'
'bean', '3803'
'beef', '9913' 
'bird', '8782'
'bird (turkey)', '9103'
'bird of prey', '3073808'
'bivalve mollusk', ' 6544'
'bos taurus taurus', '9913' 
'broiler', '9031'
'buffalo', '27592'
'bull', '9913'
'caenorhabditis elegans my316', '6239'
'calf', '9913'
'calves', '9913'
'canine', '9611'
'carrageenan', '2769'
'chacoan mara', '181543' 
'chelonoides carbonaria', '50047'
'chiloschista parishii seidenf.', '339076'
'chroicocephalus novaehollandiae (australian silver gull chick)', '2547444' 
'citellus musticus', '9996' 
'cotton rat', '42414'
'crab grass', '66017'
'crataegus sp.', '23159'
'crategus cuneata', '2741994'
'crategus monogyna', '23159'
'crow', '30420'
'ctenophthalmus teres', '185849'
'dairy herd', '9913'
'edessa sp.', '1225063'
'egret', '8899'
'epipcactis palustris', '210730'
'equus ferus caballus', '1114792'
'equus ferus caballus (thoroughbred)', '1114792' 
'eucalyptus sp.', '3932'
'felix catus domesticus', '9685'
'finch', '9133'
'fish (rainbow trout)', '8022'
'flea', '7509'
'fleas stenoponia insperata', '7509'
'fly', '7147'
'fox', '9625'
'frog', '8342'
'grasses', '4479'
'guar pulse', '3832'
'hare', '9980'
'hawk', '56259'
'indigofera argentea burm.f.', '198858'
'kalidium foliatum (pall.) moq', '224146'
'larus sp. (gull)', '8911' 
'laying hen', '9031'
'lettuce', '4236'
'lily', '4688'
'linum austriacum ssp. austriacum', '586375'
'malus domestica ''egremont russet''', '3750'
'malus domestica ''gala''', '3750'
'malus prunifolia (crab apple)', '106564'
'mammal', '40674'
'mango tree', '29780'
'marmot', '9992'
'marsupial', '9263'
'masson pine', '88730'
'melanaphis sacchari zehntner', '742174' 
'migratory bird', '8782'
'mitus arvalis', '47230'
'mole', '9373'
'mormidea sp.', '631403'
'murine', '39107'
'mus musculus c3h/orl', '10090'
'mus musculus c57bl/6j', '10090'
'mushroom', '155619'
'neopsylla setosa', '129375'
'neopsylla specialis specialis', '129375'
'nicotiana tabacum l.', '4097' 
'night heron', '8900'
'nosopsylla laeviceps', '507074'
'opossum', '38605'
'orange', '23513'
'oryza sativa cv. hwayoung', '4530'
'oryza sativa l.', '4530' 
'otter', '169417'
'oyster', '98302'
'pallasiomys meridianus', '261754' 
'pear tree', '3766'
'pelicanus rufescens', '1243782'
'pheasant', '9005'
'pigeon', '8930'
'pistacia vera l.', '55513'
'pisum sativum l.', '3888'
'porcine', '9823'
'pork', '9823'
'possum', '38609'
'poultry', '9031'
'protaetia brevitarsis seulensis larva', '438893'
'pyrus communis ''clapp''s favorite''', '23211' 
'raphanus sativus var. flamboyant 5', '3726'
'rhombomis opinus', '186474'
'rice plant', '4530'
'richardia scabra l.', '60230'
'saccharum (sugarcane)', '4546' 
'salix sp. (willow)', '40685'
'shrew', '9376'
'shrew of unidentified species', '9376'
'siskin', '1647189'
'solanum tuberosum l.', '4113'
'sophora davidii (franch.) skeels', '49839'
'sorbus aucuparia ''rowancroft coral pink''', '36599'
'sorghum bicolor (l.) moench', '4558'
'stevia rebaudiana bertoni', '55670'
'sugar cane yz08-1095', '4546'
'taurocerus sp.', '2709344'
'theobroma grandiflorum (isolate c174 resistant to witches broom)', '108881'
'trachymyrmex sp.', '34717' 
'turtle', '8459'
'vole', '337677'
'walnut', '16718'
'weever', '56735'
'wild bird', '8782'
'wild boars', '41807'
'zea mays cv. sweet belle', '4577'
'zea mays l.', '4577'
'zea mays var. rugosa (sweet corn)', '4577'

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
