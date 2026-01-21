# script for loading pre-enrichment data, modifying it, putting it into vectors
# and finally running enrichKO on it

if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  duckdb
)

# specific code for downloading MicrobiomeProfiler (pacman struggles with it)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("MicrobiomeProfiler")
library(MicrobiomeProfiler)


# connect to duckdb (much like using httr2 to connect to the NCBI API)
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

# enrichKO running
amph_enrich_result <- enrichKO(gene = amphibian_kos$ko, universe = unique_universe_kos$ko)

mammal_enrich_result <- enrichKO(gene = mammal_kos$ko, universe = unique_universe_kos$ko)

fish_enrich_result <- enrichKO(gene = fish_kos$ko, universe = unique_universe_kos$ko)

bird_enrich_result <- enrichKO(gene = bird_kos$ko, universe = unique_universe_kos$ko)

# save out to files so dont have to run again
saveRDS(amph_enrich_result, file = "02_analysis/post_annotation_analysis/enrichko_amphibian.rds")
saveRDS(mammal_enrich_result, file = "02_analysis/post_annotation_analysis/enrichko_mammal.rds")
saveRDS(fish_enrich_result, file = "02_analysis/post_annotation_analysis/enrichko_fish.rds")
saveRDS(bird_enrich_result, file = "02_analysis/post_annotation_analysis/enrichko_bird.rds")

# close connection to duckdb
dbDisconnect(con)
    