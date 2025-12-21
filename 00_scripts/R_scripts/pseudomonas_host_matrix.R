library(duckdb)
library(tidyverse)

pseudo_db <- "02_analysis/analysis_database/pseudomonas.db"

con <- dbConnect(duckdb(), here::here(pseudo_db), read_only = TRUE)
#dbListTables(con)

res <- dbSendQuery(con, "SELECT * FROM pseudomonas_annotations")
list_of_ko <- dbFetch(res)

ko_matrix <- list_of_ko %>%  count(Accession, KEGG_ko) %>%
  pivot_wider(names_from = KEGG_ko, 
              values_from = n, 
              values_fill = 0)


# host_matrix <- list_of_ko %>%  count(animal_group, KEGG_ko) %>%
#   pivot_wider(names_from = KEGG_ko, 
#               values_from = n, 
#               values_fill = 0)
# 
# 
unique_list_of_ko <- unique(list_of_ko)
# 
# unique_ko_matrix <- unique_list_of_ko %>%  count(Accession, KEGG_ko) %>%
#   pivot_wider(names_from = KEGG_ko, 
#               values_from = n, 
#               values_fill = 0)


host_matrix <- unique_list_of_ko %>%  count(animal_group, KEGG_ko) %>%
  pivot_wider(names_from = KEGG_ko, 
              values_from = n, 
              values_fill = 0) 

group_counts <- list_of_ko %>%
  distinct(Accession, animal_group) %>%
  count(animal_group, name = "n_accessions")

ko_rates <- host_matrix %>%
  pivot_longer(cols = -animal_group, 
               names_to = "KEGG_ko", 
               values_to = "count") %>%
  left_join(group_counts, by = "animal_group") %>%
  mutate(rate = count / n_accessions) %>% 
  select(-c(n_accessions, count)) %>% 
  pivot_wider(names_from = animal_group, values_from = rate)

write_delim(ko_rates, file = 'fish.cakes.csv', delim = ',')

dbDisconnect(con, shutdown = TRUE)

rm(con, res)
gc()
