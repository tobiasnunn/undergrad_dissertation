# This is a copy of the last two documents (02 and 02.1), so I am removing comments
# libraries ---------------------------------------------------------------


library(httr2)      # used to access NCBI REST API
library(tidyverse)  # all purpose tool
library(tidyjson)   # tidyverse extension for .json files
library(jsonlite)   # read in .json files
library(RPostgres)  # allows issueing of commands to postgres database
library(DBI)        # controls connection to database
library(flextable)  # tabling

# setup -------------------------------------------------------------------

# load in secrets (passwords etc) needed to link to stuff
secrets <- read.delim("../98_config/privatestuff.csv", header = T, sep = ",")
urlstring <- 'https://api.ncbi.nlm.nih.gov/datasets/v2/'
api_header <- secrets %>%
  filter(username == "ncbi") %>% pull(password) %>% as.character()
database_password <- secrets %>%
  filter(username == "postgres") %>% pull(password) %>% as.character()

new_taxa <- data.frame(taxon_id = c(561, 1763, 590, 629, 5806), genus = c("escherichia", "mycobacterium", "salmonella", "yersinia", "cryptosporidium"))

# 01. Query the NCBI API --------------------------------------------------
# get the info for escherichia down here

for (i in 1:nrow(new_taxa)) {
  # i <- 1 
  taxon <- new_taxa$taxon_id[i] 
  genus <- new_taxa$genus[i] 
  
  req <- request(urlstring) %>% 
    req_auth_bearer_token(api_header) %>%
    req_url_path_append('genome', 'taxon', taxon, 'dataset_report') %>%
    req_url_query(
      page_size = 200, 
      filters.exclude_paired_reports = TRUE) %>%
    req_headers(Accept = 'application/json') 
    
  resps <- req_perform_iterative( 
    req,  
    next_req = iterate_with_cursor( 
      'page_token', 
      resp_param = function(resp) {
        resp_body_json(resp)$next_page_token
      } 
    ),
    path = paste0("../02_analysis/02_new_genera/", genus, "_", "{Sys.Date()}_page{i}.json")  
  )
}

# 02. load output files to SQL database -----------------------------------


conn <- dbConnect(RPostgres::Postgres(), # begin connection
                  dbname = 'postgres', # name I made for database
                  host = 'localhost', # default host name (meaning my computer)
                  port = 5432, # default port for database
                  user = 'postgres', # root user 
                  password = database_password) # password

if (!dbExistsTable(conn, "genome_reports")) { 
  dbExecute(conn, paste0( 
    "CREATE TABLE ", "genome_reports", " (
        id SERIAL PRIMARY KEY,
        accession TEXT UNIQUE NOT NULL,
        genus_id TEXT NULL,
        report_data JSONB NOT NULL,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
      )"
  )) 
} 

# Initialise a progress bar
pb <- txtProgressBar(min = 0, max = nrow(new_taxa), style = 3)

# table now definitely exists, just have to populate it
for (i in 1:nrow(new_taxa)) {
  # i <- 1 # debug line
  taxon <- new_taxa$taxon_id[i] # set variables for URL
  genus <- new_taxa$genus[i] # set variable to create file names
  
  file_pattern <- paste0(genus, "_*.json") # look for files starting with spec genus name, then _(anything).json ...
  directory <- "../02_analysis/02_new_genera/" # ... in this directory
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


new_taxa_vector <- paste(new_taxa$taxon_id, collapse = ",")

req <- request(urlstring) %>% # link to NCBI database
  req_auth_bearer_token(api_header) %>% # prove I am trusted by site
  req_url_path_append('taxonomy', 'taxon', new_taxa_vector, 'dataset_report')  %>% 
  req_headers(Accept = 'application/json') %>%  # take the output in .json format
  req_perform(path = paste0("../02_analysis/02_new_genera/", Sys.Date(), "_new_taxa.json"))

dbExecute(conn, 
  "ALTER TABLE public.genus_info ADD old_new VARCHAR(10);")

dbExecute(conn, 
          "UPDATE public.genus_info SET old_new = 'existing';")

json_list <- jsonlite::read_json(paste0("../02_analysis/02_new_genera/", Sys.Date(), "_new_taxa.json"))
pb <- txtProgressBar(min = 0, max = nrow(new_taxa), style = 3)

json_reports <- json_list$reports
for (k in seq_along(json_reports)) {
  single_report <- json_reports[[k]]$taxonomy
  report_json <- jsonlite::toJSON(single_report, auto_unbox = TRUE, pretty = TRUE)
  genus_id <- single_report$tax_id
  genus_name <- single_report$current_scientific_name$name
  
  dbExecute(
    conn,
    paste0("INSERT INTO ", "genus_info", " (genus_id, genus_name, report_data, old_new) VALUES ($1, $2, $3, 'new') 
              ON CONFLICT (genus_id) DO UPDATE SET report_data = $3"), # if already exists, will update
    params = list(genus_id, genus_name, report_json))
  
  # Print progress
  setTxtProgressBar(pb, k)
}
close(pb)

# 04. summarise database data ---------------------------------------------

# so I have all the data in SQL database, and I want to transform and bring it to
# R to analyse it

# create connection to database and create table if doesnt exist already
conn <- dbConnect(RPostgres::Postgres(), # begin connection
                  dbname = 'postgres', # name I made for database
                  host = 'localhost', # default host name (meaning my computer)
                  port = 5432, # default port for database
                  user = 'postgres', # root user 
                  password = database_password) # password

result <- dbGetQuery(conn, # new command, will execute script and bring reult to "result"
                    "SELECT Q1.genus_name, Q1.host, COUNT(*) as host_count 
                    FROM 
                        (SELECT gi.genus_id, gi.genus_name, gi.old_new, gr.report_data #>> '{assembly_info, biosample, host}' as host
                        FROM public.genome_reports gr
                        INNER JOIN public.genus_info gi ON gi.genus_id = gr.genus_id) AS Q1
                    WHERE Q1.host IS NOT NULL AND
                        Q1.host NOT IN ('missing', 'mssing', 'not applicable', 'not provided', 'unknown', 'N/A', 'not collected', 'Not collected', 'no collected', 'not available', 'Not applicable')
                        AND Q1.old_new = 'new'
                    GROUP BY Q1.genus_name, Q1.host
                    ORDER BY Q1.genus_name, host_count DESC;")

# 05. turn it into a table ------------------------------------------------

table_data <- result %>%  
  group_by(genus_name) %>% # data was grouped in SQL, lost in download, needs re-grouping
  # with_ties = F means only 5 come back for each group, in the case when the values are the same for e.g. 4,5,6 and 7
  slice_max(order_by = host_count, n = 5, with_ties = F) # take top 5 (from host_content) from each group

ft <- table_data %>%
  select(genus_name, host, host_count) %>%
  flextable() %>%
  merge_v(j = "genus_name") %>%  # Merge cells vertically for genus_name
  theme_vanilla() %>%
  autofit()
ft
# Claude.ai gave me the merge_v step to combine the genus name

# save out as img
save_as_image(ft, "../04_images/new_genera_images/top_5_hosts_new.png", res = 600)


# 03 totals and proportions---------------------------------------------------------------

result <- dbGetQuery(conn,
  "SELECT gi.genus_name, gi.old_new, 'all' as count_type, COUNT(*) as total_counts FROM
  public.genome_reports gr
  INNER JOIN public.genus_info gi ON gi.genus_id = gr.genus_id
  GROUP BY gi.genus_name, gi.old_new
  UNION ALL
  SELECT gi.genus_name, gi.old_new, 'host_present' as count_type, COUNT(*) as total_counts FROM
  public.genome_reports gr
  INNER JOIN public.genus_info gi ON gi.genus_id = gr.genus_id
  WHERE gr.report_data #>> '{assembly_info, biosample, host}' IS NOT NULL AND
  gr.report_data #>> '{assembly_info, biosample, host}' NOT IN ('missing', 'mssing', 'not applicable', 'unknown', 'N/A', 'not collected', 'no collected', 'not available', 'Not applicable')
  GROUP BY gi.genus_name, gi.old_new;")
table_data <- (result %>% 
  pivot_wider(names_from = count_type, values_from = total_counts)) %>% 
  mutate(prop = host_present / all) %>% 
  arrange(desc(prop))

ft <- table_data %>%
  flextable() %>%  # Merge cells vertically for genus_name
  theme_vanilla() %>%
  bg(i = ~ old_new == "new", bg = "#fff1d9", j = 2) %>%
  bg(i = ~ old_new == "existing", bg = "#dfe7f1", j = 2) %>%
  
  autofit() %>% 
  mk_par(
    j = "prop",
    value = as_paragraph(as_chunk(prop, formatter = fmt_pct))) %>% 
  bg(j = "prop", 
     bg = scales::col_numeric(palette = "viridis", 
                              domain = c(0, 1),
                              na.color = "transparent")) %>% 
  color(~ prop < 0.5, color = "white", ~ prop) %>% 
  italic(j = "genus_name", italic = TRUE, part = "body") %>% 
  set_header_labels(genus_name = "Genus",
                    old_new = "New or Existing in dataset",
                    all = "Total accessions",
                    host_present = "Accessions with host data",
                    prop = "Proportion with host data")
ft


# save out as img
save_as_image(ft, "../04_images/new_genera_images/host_proportions_new.png", res = 600)

