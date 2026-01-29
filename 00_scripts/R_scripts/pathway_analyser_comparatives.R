# other stuff Aaron told me to do with the heatmap stuff to compare with my method (He doesnt think I have done it right)

#libraries
library(duckdb)
library(tidyverse)

# load stuff needed for enrichKO
con <- dbConnect(duckdb(), here::here("02_analysis/analysis_database/pseudomonas.db"), read_only=TRUE)

# get background universe for enrichKO (making unique explained in notebook 21)
unique_universe_kos <- dbGetQuery(con, 
                                  "SELECT DISTINCT awk.ko 
                          FROM animal_hosted_pseudo ahp
                          INNER JOIN accession_w_ko awk 
                          ON awk.accession = ahp.accession
                          WHERE ahp.animal_group IN ('Mammal', 'Fish', 'Bird', 'Amphibian');")

# amphibian kos
amphibian_kos <- dbGetQuery(con, 
                            "SELECT DISTINCT awk.ko 
                          FROM animal_hosted_pseudo ahp
                          INNER JOIN accession_w_ko awk 
                          ON awk.accession = ahp.accession
                          WHERE ahp.animal_group = 'Amphibian';")

# mammal kos
mammal_kos <- dbGetQuery(con, 
                         "SELECT DISTINCT awk.ko 
                          FROM animal_hosted_pseudo ahp
                          INNER JOIN accession_w_ko awk 
                          ON awk.accession = ahp.accession
                          WHERE ahp.animal_group = 'Mammal';")

# fish kos
fish_kos <- dbGetQuery(con, 
                       "SELECT DISTINCT awk.ko 
                          FROM animal_hosted_pseudo ahp
                          INNER JOIN accession_w_ko awk 
                          ON awk.accession = ahp.accession
                          WHERE ahp.animal_group = 'Fish';")

# bird kos
bird_kos <- dbGetQuery(con, 
                       "SELECT DISTINCT awk.ko 
                          FROM animal_hosted_pseudo ahp
                          INNER JOIN accession_w_ko awk 
                          ON awk.accession = ahp.accession
                          WHERE ahp.animal_group = 'Bird';")

# close connection to duckdb
dbDisconnect(con)

# enrichKO running
amph_enrich_result <- enrichKO(gene = amphibian_kos$ko, universe = unique_universe_kos$ko)

mammal_enrich_result <- enrichKO(gene = mammal_kos$ko, universe = unique_universe_kos$ko)

fish_enrich_result <- enrichKO(gene = fish_kos$ko, universe = unique_universe_kos$ko)

bird_enrich_result <- enrichKO(gene = bird_kos$ko, universe = unique_universe_kos$ko)

# 1. taking the whole enrichKO universe
