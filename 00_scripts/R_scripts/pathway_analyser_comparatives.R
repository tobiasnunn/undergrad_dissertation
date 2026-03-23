# other stuff Aaron told me to do with the heatmap stuff to compare with my method (He doesnt think I have done it right)

#libraries
library(duckdb)
library(tidyverse)
library(MicrobiomeProfiler)


# full universe heatmap ---------------------------------------------------


# load stuff needed for enrichKO
con <- dbConnect(duckdb(), here::here("02_analysis/analysis_database/pseudomonas.db"), read_only=TRUE)

# amphibian kos
amphibian_kos <- dbGetQuery(con, 
                            "SELECT DISTINCT pa.KEGG_ko 
                          FROM pseudomonas_annotations pa
                          WHERE pa.animal_group = 'Amphibian';")

# mammal kos
mammal_kos <- dbGetQuery(con, 
                         "SELECT DISTINCT pa.KEGG_ko 
                          FROM pseudomonas_annotations pa
                          WHERE pa.animal_group = 'Mammal';")

# fish kos
fish_kos <- dbGetQuery(con, 
                       "SELECT DISTINCT pa.KEGG_ko 
                          FROM pseudomonas_annotations pa
                          WHERE pa.animal_group = 'Fish';")

# bird kos
bird_kos <- dbGetQuery(con, 
                       "SELECT DISTINCT pa.KEGG_ko 
                          FROM pseudomonas_annotations pa
                          WHERE pa.animal_group = 'Bird';")

# close connection to duckdb
dbDisconnect(con)

# 1. taking the whole enrichKO universe

# enrichKO running
amph_enrich_no_uni <- enrichKO(gene = amphibian_kos$KEGG_ko)

mammal_enrich_no_uni <- enrichKO(gene = mammal_kos$KEGG_ko)

fish_enrich_no_uni <- enrichKO(gene = fish_kos$KEGG_ko)

bird_enrich_no_uni <- enrichKO(gene = bird_kos$KEGG_ko)

# extract result object containing test information
amphstats <- data.frame(amph_enrich_no_uni@result)
mamstats <- data.frame(mammal_enrich_no_uni@result)
fishstats <- data.frame(fish_enrich_no_uni@result)
birdstats <- data.frame(bird_enrich_no_uni@result)

# combine objects to allow filtering for comparisons
all_groups <- amphstats %>% 
  full_join(mamstats, by = c("ID", "Description", "BgRatio"), 
            suffix = c(".amphibian", ".mammal")) %>%
  full_join(birdstats, by = c("ID", "Description", "BgRatio")) %>%
  full_join(fishstats, by = c("ID", "Description", "BgRatio"), suffix = c(".bird", ".fish"))

# find the differentially expressed pathways for each group (probably a cleaner way to do it)
#amphibians
differential_pathways_amph <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                       all_groups$p.adjust.fish > 0.05 &
                                       all_groups$p.adjust.mammal > 0.05 &
                                       all_groups$p.adjust.bird > 0.05)
#none

#mammals
differential_pathways_mam <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                      all_groups$p.adjust.fish > 0.05 &
                                      all_groups$p.adjust.mammal < 0.05 &
                                      all_groups$p.adjust.bird > 0.05)
# 4 - mostly about metabolism, no share with ones found in previous analysis

#birds
differential_pathways_bird <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                       all_groups$p.adjust.fish > 0.05 &
                                       all_groups$p.adjust.mammal > 0.05 &
                                       all_groups$p.adjust.bird < 0.05)
# none

#fish
differential_pathways_fish <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                       all_groups$p.adjust.fish < 0.05 &
                                       all_groups$p.adjust.mammal > 0.05 &
                                       all_groups$p.adjust.bird > 0.05)
# 2 - more metabolism, no overlap with previous analysis

# so... still got six, but a different six and now the groups have flipped, so that instead of amph and bird, now the pathways come from mam and fish. what the actual hell?


#combine them
differential_pathways <- rbind(differential_pathways_amph, differential_pathways_mam, differential_pathways_fish, differential_pathways_bird)

differential_pathways

# create the base heatmap source item
heatmap_source <- differential_pathways %>%  select(c(ID, p.adjust.amphibian,
                                           p.adjust.mammal, p.adjust.bird,
                                           p.adjust.fish))
colnames(heatmap_source) <- c("ID", "Amphibian", "Mammal", "Bird", "Fish")
# make it long
heat_source_long <- heatmap_source %>% 
  pivot_longer(cols = c(Amphibian,
                        Mammal, Bird,
                        Fish), names_to = "p.adjusted", 
               values_to = "p_value")

# add colour to do that 1 square in white text
heat_source_long <- heat_source_long %>%
  mutate(text_color = ifelse(p_value > 0.5, "white", "black"))

# basic plot
base <- ggplot(heat_source_long, aes(ID, p.adjusted, fill= p_value)) + 
  geom_tile() +
  scale_fill_viridis_c(direction = -1) +
  theme_minimal()

# add border to significant cells (where p < 0.05)
# got this off of google ai, realise I probably shouldnt have, but it works
developed <- base +
  # Add border to cells > threshold
  geom_tile(data = filter(heat_source_long, p_value < 0.05), 
            color = "black", size = 1, fill = NA) +
  labs(x = "KEGG Pathway ID", y = "Host Taxon", fill = "Adjusted\nP value") +
  geom_text(data = heat_source_long,
            aes(label = round(p_value, digits = 3), color = text_color), size = 3) +
  scale_color_identity()
developed
# save
ggsave(filename = "03_outputs/02_enrichment/full_universe_heatmap.jpeg", plot = developed, width = 20, height = 15, units = "cm")


# DEVO - REDUNDANT --------------------------------------------------------


#------------# 
# shared universe pathways #
# pathways where all 4 are enriched
differential_pathways_all <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                      all_groups$p.adjust.fish < 0.05 &
                                      all_groups$p.adjust.mammal < 0.05 &
                                      all_groups$p.adjust.bird < 0.05)


# pathways where 3 are enriched and 1 isnt
differential_pathways_noamph <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                         all_groups$p.adjust.fish < 0.05 &
                                         all_groups$p.adjust.mammal < 0.05 &
                                         all_groups$p.adjust.bird < 0.05) 

differential_pathways_nofish <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                         all_groups$p.adjust.fish > 0.05 &
                                         all_groups$p.adjust.mammal < 0.05 &
                                         all_groups$p.adjust.bird < 0.05) 

differential_pathways_nomam <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                        all_groups$p.adjust.fish < 0.05 &
                                        all_groups$p.adjust.mammal > 0.05 &
                                        all_groups$p.adjust.bird < 0.05) 

differential_pathways_nobird <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                         all_groups$p.adjust.fish < 0.05 &
                                         all_groups$p.adjust.mammal < 0.05 &
                                         all_groups$p.adjust.bird > 0.05) 
# no pathways where 3 are enriched and 1 isnt

# 2 on 2 off
# amphibian + other
differential_pathways_amphfish <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                           all_groups$p.adjust.fish < 0.05 &
                                           all_groups$p.adjust.mammal > 0.05 &
                                           all_groups$p.adjust.bird > 0.05) 

differential_pathways_amphmam <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                          all_groups$p.adjust.fish > 0.05 &
                                          all_groups$p.adjust.mammal < 0.05 &
                                          all_groups$p.adjust.bird > 0.05)

differential_pathways_amphbird <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                           all_groups$p.adjust.fish > 0.05 &
                                           all_groups$p.adjust.mammal > 0.05 &
                                           all_groups$p.adjust.bird < 0.05)
# fish + other

# amphibian already done

differential_pathways_fishmam <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                          all_groups$p.adjust.fish < 0.05 &
                                          all_groups$p.adjust.mammal < 0.05 &
                                          all_groups$p.adjust.bird > 0.05)

differential_pathways_fishbird <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                           all_groups$p.adjust.fish < 0.05 &
                                           all_groups$p.adjust.mammal > 0.05 &
                                           all_groups$p.adjust.bird < 0.05)
# mammal + other

# amphibian already done
# fish already done

differential_pathways_mambird <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                          all_groups$p.adjust.fish > 0.05 &
                                          all_groups$p.adjust.mammal < 0.05 &
                                          all_groups$p.adjust.bird < 0.05)

# bird + other
# amph done / fish done / mam done


#combine them
shared_differential_pathways <- rbind(differential_pathways_all,
                                      differential_pathways_amphbird, 
                                      differential_pathways_amphfish,
                                      differential_pathways_amphmam,
                                      differential_pathways_fishbird,
                                      differential_pathways_fishmam,
                                      differential_pathways_mambird,
                                      differential_pathways_noamph,
                                      differential_pathways_nobird,
                                      differential_pathways_nofish,
                                      differential_pathways_nomam)

heatmap_source <- shared_differential_pathways %>%  select(c(ID, p.adjust.amphibian,
                                                      p.adjust.mammal, p.adjust.bird,
                                                      p.adjust.fish))
colnames(heatmap_source) <- c("ID", "Amphibian", "Mammal", "Bird", "Fish")
# make it long
heat_source_long <- heatmap_source %>% 
  pivot_longer(cols = c(Amphibian,
                        Mammal, Bird,
                        Fish), names_to = "p.adjusted", 
               values_to = "p_value")

# add colour to do that 1 square in white text
heat_source_long <- heat_source_long %>%
  mutate(text_color = ifelse(p_value > 0.1, "white", "black"))

# basic plot
base <- ggplot(heat_source_long, aes(ID, p.adjusted, fill= p_value)) + 
  geom_tile() +
  scale_fill_viridis_c(direction = -1) +
  theme_minimal()

# add border to significant cells (where p < 0.05)
# got this off of google ai, realise I probably shouldnt have, but it works
developed <- base +
  # Add border to cells > threshold
  geom_tile(data = filter(heat_source_long, p_value < 0.05), 
            color = "black", size = 1, fill = NA) +
  labs(x = "KEGG Pathway ID", y = "Host Taxon", fill = "Adjusted\nP value") +
  geom_text(data = heat_source_long,
            aes(label = round(p_value, digits = 3), color = text_color), size = 3) +
  scale_color_identity() + coord_flip()
developed
# save
ggsave(filename = "03_outputs/02_enrichment/shared_universe_heatmap.jpeg", plot = developed, width = 20, height = 40, units = "cm")
# pathways where some / all have it ---------------------------------------

# from "enrichko_runner.R"

# amph_paths <- readRDS("02_analysis/post_annotation_analysis/enrichko_amphibian.rds")
# mam_paths <- readRDS("02_analysis/post_annotation_analysis/enrichko_mammal.rds")
# fish_paths <- readRDS("02_analysis/post_annotation_analysis/enrichko_fish.rds")
# bird_paths <- readRDS("02_analysis/post_annotation_analysis/enrichko_bird.rds")
con <- dbConnect(duckdb(), here::here("02_analysis/analysis_database/pseudomonas.db"), read_only=TRUE)

# get background universe for enrichKO (making unique explained in notebook 21)
unique_universe_kos <- dbGetQuery(con, 
                                  "SELECT * 
                          FROM pseudomonas_annotations pa;")

# amphibian kos
amphibian_kos <- dbGetQuery(con, 
                            "SELECT DISTINCT pa.KEGG_ko 
                          FROM pseudomonas_annotations pa
                          WHERE pa.animal_group = 'Amphibian';")

# mammal kos
mammal_kos <- dbGetQuery(con, 
                         "SELECT DISTINCT pa.KEGG_ko 
                          FROM pseudomonas_annotations pa
                          WHERE pa.animal_group = 'Mammal';")

# fish kos
fish_kos <- dbGetQuery(con, 
                       "SELECT DISTINCT pa.KEGG_ko 
                          FROM pseudomonas_annotations pa
                          WHERE pa.animal_group = 'Fish';")

# bird kos
bird_kos <- dbGetQuery(con, 
                       "SELECT DISTINCT pa.KEGG_ko 
                          FROM pseudomonas_annotations pa
                          WHERE pa.animal_group = 'Bird';")

# enrichKO running
amph_enrich_result <- enrichKO(gene = amphibian_kos$KEGG_ko, universe = unique_universe_kos$KEGG_ko)

mammal_enrich_result <- enrichKO(gene = mammal_kos$KEGG_ko, universe = unique_universe_kos$KEGG_ko)

fish_enrich_result <- enrichKO(gene = fish_kos$KEGG_ko, universe = unique_universe_kos$KEGG_ko)

bird_enrich_result <- enrichKO(gene = bird_kos$KEGG_ko, universe = unique_universe_kos$KEGG_ko)

dbDisconnect(con)

# extract result object containing test information
amphstats <- data.frame(amph_enrich_result@result)
mamstats <- data.frame(mammal_enrich_result@result)
fishstats <- data.frame(fish_enrich_result@result)
birdstats <- data.frame(bird_enrich_result@result)

# combine objects to allow filtering for comparisons
all_groups <- amphstats %>% #as description and bgratio relate to the pathway and not the group, they are de-duplicated by joining on them
  full_join(mamstats, by = c("ID", "Description", "BgRatio"), suffix = c(".amphibian", ".mammal")) %>%
  full_join(birdstats, by = c("ID", "Description", "BgRatio")) %>% # trying to add suffix does nothing as there is no conflict (notebook 22)
  full_join(fishstats, by = c("ID", "Description", "BgRatio"), suffix = c(".bird", ".fish"))

# pathways where all 4 are enriched
differential_pathways_all <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                       all_groups$p.adjust.fish < 0.05 &
                                       all_groups$p.adjust.mammal < 0.05 &
                                       all_groups$p.adjust.bird < 0.05)


# pathways where 3 are enriched and 1 isnt
differential_pathways_noamph <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                      all_groups$p.adjust.fish < 0.05 &
                                      all_groups$p.adjust.mammal < 0.05 &
                                      all_groups$p.adjust.bird < 0.05) 

differential_pathways_nofish <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                         all_groups$p.adjust.fish > 0.05 &
                                         all_groups$p.adjust.mammal < 0.05 &
                                         all_groups$p.adjust.bird < 0.05) 

differential_pathways_nomam <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                         all_groups$p.adjust.fish < 0.05 &
                                         all_groups$p.adjust.mammal > 0.05 &
                                         all_groups$p.adjust.bird < 0.05) 

differential_pathways_nobird <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                         all_groups$p.adjust.fish < 0.05 &
                                         all_groups$p.adjust.mammal < 0.05 &
                                         all_groups$p.adjust.bird > 0.05) 
# no pathways where 3 are enriched and 1 isnt

# 2 on 2 off
# amphibian + other
differential_pathways_amphfish <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                      all_groups$p.adjust.fish < 0.05 &
                                      all_groups$p.adjust.mammal > 0.05 &
                                      all_groups$p.adjust.bird > 0.05) 

differential_pathways_amphmam <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                           all_groups$p.adjust.fish > 0.05 &
                                           all_groups$p.adjust.mammal < 0.05 &
                                           all_groups$p.adjust.bird > 0.05)

differential_pathways_amphbird <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                          all_groups$p.adjust.fish > 0.05 &
                                          all_groups$p.adjust.mammal > 0.05 &
                                          all_groups$p.adjust.bird < 0.05)
# fish + other

# amphibian already done

differential_pathways_fishmam <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                           all_groups$p.adjust.fish < 0.05 &
                                           all_groups$p.adjust.mammal < 0.05 &
                                           all_groups$p.adjust.bird > 0.05)

differential_pathways_fishbird <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                          all_groups$p.adjust.fish < 0.05 &
                                          all_groups$p.adjust.mammal > 0.05 &
                                          all_groups$p.adjust.bird < 0.05)
# mammal + other

# amphibian already done
# fish already done

differential_pathways_mambird <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                          all_groups$p.adjust.fish > 0.05 &
                                          all_groups$p.adjust.mammal < 0.05 &
                                          all_groups$p.adjust.bird < 0.05)

# bird + other
# amph done / fish done / mam done


#combine them
shared_differential_pathways <- rbind(differential_pathways_all,
                                      differential_pathways_amphbird, 
                                      differential_pathways_amphfish,
                                      differential_pathways_amphmam,
                                      differential_pathways_fishbird,
                                      differential_pathways_fishmam,
                                      differential_pathways_mambird,
                                      differential_pathways_noamph,
                                      differential_pathways_nobird,
                                      differential_pathways_nofish,
                                      differential_pathways_nomam)

# create the base heatmap source item
heatmap_source <- shared_differential_pathways %>%  select(c(ID, p.adjust.amphibian,
                                                      p.adjust.mammal, p.adjust.bird,
                                                      p.adjust.fish))
colnames(heatmap_source) <- c("ID", "Amphibian", "Mammal", "Bird", "Fish")

# make it long
heat_source_long <- heatmap_source %>% 
  pivot_longer(cols = c(Amphibian,
                        Mammal, Bird,
                        Fish), names_to = "p.adjusted", 
               values_to = "p_value")

# add colour to do that 1 square in white text
heat_source_long <- heat_source_long %>%
  mutate(text_color = ifelse(p_value > 0.2, "white", "black"))

# basic plot
base <- ggplot(heat_source_long, aes(ID, p.adjusted, fill= p_value)) + 
  geom_tile() +
  scale_fill_viridis_c(direction = -1) +
  theme_minimal()

# add border to significant cells (where p < 0.05)
# got this off of google ai, realise I probably shouldnt have, but it works
developed <- base +
  # Add border to cells > threshold
  geom_tile(data = filter(heat_source_long, p_value < 0.05), 
            color = "black", size = 1, fill = NA) +
  labs(x = "KEGG Pathway ID", y = "Host Taxon", fill = "Adjusted\nP value") +
  geom_text(data = heat_source_long,
            aes(label = round(p_value, digits = 3), color = text_color), size = 3) +
  scale_color_identity()
developed
# save
ggsave(filename = "03_outputs/02_enrichment/shared_differences_heatmap.jpeg", plot = developed, width = 15, height = 20, units = "cm")



# Aarons method -----------------------------------------------------------

# enrich individual by individual and take significant correlation in a group
# as where 1 group has over 80% of individuals significantly enriched, and
# the others sum to less than 50%

# load stuff needed for enrichKO
con <- dbConnect(duckdb(), here::here("02_analysis/analysis_database/pseudomonas.db"), read_only=TRUE)

# bring in whole object
pseudomonas_annotations <- dbGetQuery(con, 
                            "SELECT * 
                          FROM pseudomonas_annotations pa;")



# close connection to duckdb
dbDisconnect(con)

# get unique list of accessions to use in loop
accession_list <- unique(pseudomonas_annotations$Accession)

# 1. run enrichKO in a for loop so it enriches per individual

# create hollow list
enriched_individuals <- list()
# set i for naming values in list
i <- 1

# enrichKO running
for(accession in accession_list){
  # debug: accession <- "GCA_025792185.1"
  # check I can do more than 1: accession <- "GCA_037967095.1"
  # name the object for the list
  name <- accession
  # get the value from enrichKO
  value <- enrichKO(gene = pseudomonas_annotations[pseudomonas_annotations$Accession == accession, "KEGG_ko"])
  # add the value to the list
  enriched_individuals <- append(enriched_individuals, values = value)
  # modify the name so I know what value relates to what accession
  names(enriched_individuals)[i] <- name
  i <- i + 1
}

# write object out to RDS so I dont have to rerun
saveRDS(enriched_individuals, file = "02_analysis/post_annotation_analysis/enrich_by_individual.rds")

enriched_individuals <- readRDS("02_analysis/post_annotation_analysis/enrich_by_individual.rds")
### create a dataframe of the key information ###

# 2 parts, the host info and the enrich info, join on accession?
# host info
pseudo_info <- unique(data.frame("Host_group" = pseudomonas_annotations$animal_group,
                          "Accession" = pseudomonas_annotations$Accession))

# enrich info
y <- data.frame()

for (individual in names(enriched_individuals)) {
  z <- enriched_individuals[[individual]]@result
  z$Accession <- individual
  y <- rbind(y, z)
}

# write it out
saveRDS(y, file = "02_analysis/post_annotation_analysis/enriched_by_accession.rds")

# read out to duckdb
con <- dbConnect(duckdb(), 
                 dbdir = here::here("02_analysis/analysis_database/pseudomonas.db"),
                 read_only = F)

duckdb::dbWriteTable(con, "enriched_by_accession", y)

dbDisconnect(con)

# did some stuff in SQL, but the code is going to be used here to call it and
# get the table from duckdb

saveRDS(y, file = "02_analysis/post_annotation_analysis/enriched_by_accession.rds")

# read out to duckdb
con <- dbConnect(duckdb(), 
                 dbdir = here::here("02_analysis/analysis_database/pseudomonas.db"),
                 read_only = F)

enrichment_by_indv <- dbGetQuery(con, 'WITH distinct_animal_group AS 
(SELECT DISTINCT Accession, animal_group FROM
  pseudomonas_annotations),
animal_group_count AS
(SELECT animal_group, COUNT(*) AS total_accessions
  FROM distinct_animal_group 
  GROUP BY animal_group ),
significant_pathways AS 
(SELECT eba.Accession, eba.ID, dag.animal_group 
  FROM enriched_by_accession eba
  LEFT JOIN distinct_animal_group dag
  ON dag.Accession = eba.Accession
  WHERE "p.adjust" <= 0.05),
raw_proportions AS 
(SELECT sp.animal_group, sp.ID, agc.total_accessions, COUNT(*) AS enriched_accessions,
  enriched_accessions/agc.total_accessions AS proportion 
  FROM significant_pathways sp
  LEFT JOIN animal_group_count agc
  ON agc.animal_group = sp.animal_group 
  GROUP BY sp.animal_group, sp.ID, agc.total_accessions),
pre_pivot_proportions AS
(SELECT animal_group, ID, proportion
  FROM raw_proportions),
pivot_proportions AS 
(PIVOT pre_pivot_proportions
  ON animal_group
  USING SUM(proportion)
  GROUP BY ID)
SELECT *
FROM pivot_proportions;')


dbDisconnect(con)

# code used to identify pathways:
# SELECT ID, SUM(CASE WHEN proportion > 0.8 THEN 1 ELSE 0 END) AS high_count,
# SUM(CASE WHEN proportion <= 0.5 THEN 1 ELSE 0 END) AS low_count
# FROM pre_pivot_proportions
# GROUP BY ID
# HAVING high_count = 1 AND low_count = 3;

# -- 2 pathways out of 107 for the 1 over 80 rest under 50 analysis, not a good sign of adaptation map00930 and map00405
# -- 69 pathways out of 107 for the all above 80, lending iteself to the conclusion that 
# -- there is no real adaptation, it is all generalism, how would I visualise, or communicate that in results?

# heatmapping

# make it long
heat_source_long <- enrichment_by_indv %>% 
  pivot_longer(cols = c(Amphibian,
                        Mammal, Bird,
                        Fish), names_to = "host_group", 
               values_to = "proportion")

# turn NAs (how some part of the process handled 0s to 0s)
heat_source_long[is.na(heat_source_long)] <- 0

# add colour to do that 1 square in white text
heat_source_long <- heat_source_long %>%
  mutate(text_color = ifelse(proportion < 0.4, "white", "black"))

# basic plot
base <- ggplot(heat_source_long, aes(ID, host_group, fill= proportion)) + 
  geom_tile() +
  scale_fill_viridis_c(direction = 1) +
  theme_minimal()

# add border to significant cells (where p < 0.05)
# got this off of google ai, realise I probably shouldnt have, but it works
developed <- base +
  # Add border to cells > threshold
  geom_tile(data = filter(heat_source_long, proportion > 0.8), 
            color = "black", size = 1, fill = NA) +
  labs(x = "KEGG Pathway ID", y = "Host Taxon", fill = "Proportion of\nindividuals") +
  geom_text(data = heat_source_long,
            aes(label = round(proportion, digits = 3), color = text_color), size = 3) +
  scale_color_identity() + coord_flip()
developed
# save
ggsave(filename = "03_outputs/02_enrichment/aaron_method_full.jpeg", plot = developed, width = 15, height = 50, units = "cm")

#------------------------------#
# just the 2 important pathways

# turn NAs (how some part of the process handled 0s to 0s)
heat_source_long[is.na(heat_source_long)] <- 0

# make it just the two pathways
heat_source_long <- heat_source_long[heat_source_long$ID %in% c("map00930", "map00405"),]

# add colour to do that 1 square in white text
heat_source_long <- heat_source_long %>%
  mutate(text_color = ifelse(proportion < 0.4, "white", "black"))

# basic plot
base <- ggplot(heat_source_long, aes(ID, host_group, fill= proportion)) + 
  geom_tile() +
  scale_fill_viridis_c(direction = 1) +
  theme_minimal()

# add border to significant cells (where p < 0.05)
# got this off of google ai, realise I probably shouldnt have, but it works
developed <- base +
  # Add border to cells > threshold
  geom_tile(data = filter(heat_source_long, proportion > 0.8), 
            color = "black", size = 1, fill = NA) +
  labs(x = "KEGG Pathway ID", y = "Host Taxon", fill = "Proportion of\nindividuals") +
  geom_text(data = heat_source_long,
            aes(label = round(proportion, digits = 3), color = text_color), size = 3) +
  scale_color_identity()
developed
# save
ggsave(filename = "03_outputs/02_enrichment/aaron_method_diffs.jpeg", plot = developed, width = 15, height = 15, units = "cm")
