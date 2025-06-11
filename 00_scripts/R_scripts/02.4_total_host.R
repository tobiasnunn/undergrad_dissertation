# This is a copy of the last two documents (02 and 02.1), so I am removing comments
# libraries ---------------------------------------------------------------


library(tidyverse)  # all purpose tool
library(RPostgres)  # allows issueing of commands to postgres database
library(DBI)        # controls connection to database
library(gt)         # more flexible tabling

# setup -------------------------------------------------------------------

# load in secrets (passwords etc) needed to link to stuff
secrets <- read.delim("../../98_config/privatestuff.csv", header = T, sep = ",")
# path needs to be this way so that the file can be run in the notebook
# subsequently it won't run on its own

database_password <- secrets %>%
  filter(username == "postgres") %>% pull(password) %>% as.character()

# 01. paginated data ------------------------------------------------------

conn <- dbConnect(RPostgres::Postgres(), 
                  dbname = 'postgres',
                  host = 'localhost', 
                  port = 5432, 
                  user = 'postgres', 
                  password = database_password) 

# get data from database
result <- dbGetQuery(conn, "SELECT Q1.host, COUNT(*) as host_count FROM 
(SELECT gi.genus_id, gi.genus_name, gr.report_data #>> '{assembly_info, biosample, host}' as host
  FROM public.genome_reports gr
  INNER JOIN public.genus_info gi ON gi.genus_id = gr.genus_id) AS Q1
WHERE Q1.host IS NOT NULL AND
Q1.host NOT IN ('missing', 'mssing', 'not applicable', 'unknown', 'N/A', 'not collected', 'no collected', 'not available', 'Not applicable')
GROUP BY Q1.host
ORDER BY host_count DESC;")

#table is generated in notebook