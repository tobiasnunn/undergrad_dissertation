# other stuff Aaron told me to do with the heatmap stuff to compare with my method (He doesnt think I have done it right)

#libraries
library(duckdb)
library(tidyverse)
library(MicrobiomeProfiler)


# full universe heatmap ---------------------------------------------------


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

# 1. taking the whole enrichKO universe

# enrichKO running
amph_enrich_no_uni <- enrichKO(gene = amphibian_kos$ko)

mammal_enrich_no_uni <- enrichKO(gene = mammal_kos$ko)

fish_enrich_no_uni <- enrichKO(gene = fish_kos$ko)

bird_enrich_no_uni <- enrichKO(gene = bird_kos$ko)

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
  mutate(text_color = ifelse(p_value > 0.8, "white", "black"))

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


# pathways where some / all have it ---------------------------------------

# from "enrichko_runner.R"

amph_paths <- readRDS("02_analysis/post_annotation_analysis/enrichko_amphibian.rds")
mam_paths <- readRDS("02_analysis/post_annotation_analysis/enrichko_mammal.rds")
fish_paths <- readRDS("02_analysis/post_annotation_analysis/enrichko_fish.rds")
bird_paths <- readRDS("02_analysis/post_annotation_analysis/enrichko_bird.rds")

# extract result object containing test information
amphstats <- data.frame(amph_paths@result)
mamstats <- data.frame(mam_paths@result)
fishstats <- data.frame(fish_paths@result)
birdstats <- data.frame(bird_paths@result)

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
shared_differential_pathways <- rbind(differential_pathways_amphbird, 
                                      differential_pathways_mambird,
                                      differential_pathways_noamph,
                                      differential_pathways_nofish)

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
ggsave(filename = "03_outputs/02_enrichment/shared_differences_heatmap.jpeg", plot = developed, width = 20, height = 15, units = "cm")
