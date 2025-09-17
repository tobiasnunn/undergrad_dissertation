
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
database_password <- secrets %>%
  filter(username == "postgres") %>% pull(password) %>% as.character()

# 01. summarise database data ---------------------------------------------

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
                        (SELECT gi.genus_id, gi.genus_name, gr.report_data #>> '{assembly_info, biosample, host}' as host
                        FROM public.genome_reports gr
                        INNER JOIN public.genus_info gi ON gi.genus_id = gr.genus_id) AS Q1
                    WHERE Q1.host IS NOT NULL AND
                        Q1.host NOT IN ('missing', 'mssing', 'not applicable', 'not provided', 'unknown', 'N/A', 'not collected', 'Not collected', 'no collected', 'not available', 'Not applicable')
                    GROUP BY Q1.genus_name, Q1.host
                    ORDER BY Q1.genus_name, host_count DESC;")
# bit of a complex script, major things that might not make sense at first:
# one can feed one query into another by using "AS Q1" - this means one does not have to use the long table paths, a "subquery", specifically
# I would have to write "gr.report_data #>> '{assembly_info, biosample, host}'", any time I wanted to refer to host
# gi and gr are ways of setting shorthands for tables, also to shorten paths for simplicity, they get set on lines 36 and 37
# the whole query summarises the number of hosts by genus, giving me a count for each value in the column, ordered by host
# I did get some SQL help for this, as one does not get this good at SQL immediately

# 02. turn it into a table ------------------------------------------------

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
save_as_image(ft, "../04_images/existing_genera_images/top_5_hosts.png", res = 600)


# 03 totals and proportions---------------------------------------------------------------

result <- dbGetQuery(conn,
  "SELECT gi.genus_name, 'all' as count_type, COUNT(*) as total_counts FROM
  public.genome_reports gr
  INNER JOIN public.genus_info gi ON gi.genus_id = gr.genus_id
  GROUP BY gi.genus_name
  UNION ALL
  SELECT gi.genus_name, 'host_present' as count_type, COUNT(*) as total_counts FROM
  public.genome_reports gr
  INNER JOIN public.genus_info gi ON gi.genus_id = gr.genus_id
  WHERE gr.report_data #>> '{assembly_info, biosample, host}' IS NOT NULL AND
  gr.report_data #>> '{assembly_info, biosample, host}' NOT IN ('missing', 'mssing', 'not applicable', 'unknown', 'N/A', 'not collected', 'no collected', 'not available', 'Not applicable')
  GROUP BY gi.genus_name;")
table_data <- (result %>% 
  pivot_wider(names_from = count_type, values_from = total_counts)) %>% 
  mutate(prop = host_present / all)

ft <- table_data %>%
  flextable() %>%  # Merge cells vertically for genus_name
  theme_vanilla() %>%
  autofit() %>% 
  mk_par(
    j = "prop",
    value = as_paragraph(as_chunk(prop, formatter = fmt_pct)))
ft

# save out as img
save_as_image(ft, "../04_images/existing_genera_images/host_proportions.png", res = 600)

