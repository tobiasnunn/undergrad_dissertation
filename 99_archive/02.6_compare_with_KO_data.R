
# libraries ---------------------------------------------------------------

library(tidyverse)  # all purpose tool
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

conn <- dbConnect(RPostgres::Postgres(), # begin connection
                  dbname = 'postgres', # name I made for database
                  host = 'localhost', # default host name (meaning my computer)
                  port = 5432, # default port for database
                  user = 'postgres', # root user 
                  password = database_password) # password


# 01. read accessions for existing KO -------------------------------------

KO_files <- data.frame(accession = list.files(path = "../01_inputs/02_previous_analysis_results/eggnog-mapper/", pattern = "*ko_pathway.tsv"))
KO_files$accession <- gsub("_ko_pathway.tsv", "", KO_files$accession)
KO_files$accession_number <- gsub("GCA_|GCF_", "", KO_files$accession)


accession_with_host <- read_delim("../03_outputs/host_evaluation/accession_to_host_class.tsv")

accession_with_host <- accession_with_host %>% 
  mutate(accession_number = gsub("GCA_|GCF_", "", accession))


accession_with_host_KO <- accession_with_host %>% 
  filter(accession %in% KO_files$accession)

accession_with_host_KO_just_number <- accession_with_host %>% 
  filter(accession_number %in% KO_files$accession_number)

#display
accession_with_host_KO_just_number_for_display <- accession_with_host_KO_just_number %>% 
  group_by(class_name, genus_name) %>% 
  count(name = "count_of_accessions") %>% 
  arrange(desc(count_of_accessions))

# now pivoted
pivoted_accession_With_KO_for_display <- accession_with_host_KO_just_number_for_display %>% 
  pivot_wider(names_from = genus_name, values_from = count_of_accessions)

write_delim(accession_with_host_KO_just_number, "../03_outputs/host_evaluation/accession_to_class_KO.tsv", delim = "\t")

# should have requested the primary identifier
# there are only 5 instances where the ones I have used is different
# TODO: take the GCF/A off the front in order to do the quick-fix
# TODO: update minor issues to include the inconsistency in accession GCF/A and it might 
# be hard to get a clean way to do that, might require some redo-ing, aybe use an alias table
# TODO: get this info to the notebook in one way or another