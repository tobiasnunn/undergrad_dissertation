
# libraries ---------------------------------------------------------------
library(dplyr)

# bring in enrichKO output rds files --------------------------------------
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
all_groups <- amphstats %>% #as description and bgratio relate to the pathway and hot the group, they are de-duplicated by joining on them
  full_join(mamstats, by = c("ID", "Description", "BgRatio"), suffix = c(".amphibian", ".mammal")) %>%
  full_join(birdstats, by = c("ID", "Description", "BgRatio")) %>% # trying to add suffix does nothing as there is no conflict (notebook 22)
  full_join(fishstats, by = c("ID", "Description", "BgRatio"), suffix = c(".bird", ".fish"))

# find the differentially expressed pathways for each group (probably a cleaner way to do it)
#amphibians
differential_pathways_amph <- filter(all_groups, all_groups$p.adjust.amphibian < 0.05 & 
                                  all_groups$p.adjust.fish > 0.05 &
                                  all_groups$p.adjust.mammal > 0.05 &
                                  all_groups$p.adjust.bird > 0.05) # this is birds, need to fix that step so it has the correct suffix

#mammals
differential_pathways_mam <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                       all_groups$p.adjust.fish > 0.05 &
                                       all_groups$p.adjust.mammal < 0.05 &
                                       all_groups$p.adjust.bird > 0.05) # this is birds, need to fix that step so it has the correct suffix
#hmmm no mammal pathways, notes of lack of specialisation, mammals too broad? Expected as make up 91% of background. perhaps I should do an inside-mammal analysis, instead of putting it as a "next step" of the research.

#birds
differential_pathways_bird <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                       all_groups$p.adjust.fish > 0.05 &
                                       all_groups$p.adjust.mammal > 0.05 &
                                       all_groups$p.adjust.bird < 0.05) # this is birds, need to fix that step so it has the correct suffix

#fish
differential_pathways_fish <- filter(all_groups, all_groups$p.adjust.amphibian > 0.05 & 
                                       all_groups$p.adjust.fish < 0.05 &
                                       all_groups$p.adjust.mammal > 0.05 &
                                       all_groups$p.adjust.bird > 0.05) # this is birds, need to fix that step so it has the correct suffix
# hmm, none from fish either, I remember when I did a little look into the geography of the hosts, that all the fish come from 1 sampling event in like a river in China, might explain it, maybe because not broad enough range of hosts, pathways that are enriched in fish have been missed by chance, maybe analysis is unfair, as it is too specific (i.e. mammal vs amphibian vs bird vs "chinies river salmon")

# Where did I do that little geography test, and was it part of that 80% vs 50% thing I did for a bit of fun a while back too?

#combine them
differential_pathways <- rbind(differential_pathways_amph, differential_pathways_mam, differential_pathways_fish, differential_pathways_bird)

# write to RDS (so I dont have to run all this code again once I clear the environment to work on the phylogenetic tree)
saveRDS(differential_pathways, file = "02_analysis/post_annotation_analysis/differential_pathways.rds")

# read that back in when making a heatmap in `heatmap_maker.R`


# testing the pre-suffixing method of suffixing ---------------------------

# # testing the pre-suffixing method of suffixing
# z_amphstats <- data.frame(amph_paths@result)
# #got stuck and asked Claude as I dont have time to spend hours mucking about with this
# # base R
# cols_to_rename <- !colnames(z_amphstats) %in% c("ID", "Description")
# colnames(z_amphstats)[cols_to_rename] <- paste0(colnames(z_amphstats)[cols_to_rename], ".amph")
# #Tidyverse version (shorter but requires understanding ~):
# rz_amphstats <- z_amphstats %>% 
#   rename_with(~paste0(., ".amph"), -c(ID, Description))
# colnames(z_amphstats[,c(3,5:12)]) <- c(z_amphstats[,c(3,5:12)], ".amph")
