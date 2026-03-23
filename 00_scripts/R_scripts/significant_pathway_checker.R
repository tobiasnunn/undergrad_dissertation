# explanation -------------------------------------------------------------
# I want to check my results to see if they are significant the way I did it
# Specifically, I want to see if the pathways enriched are ghost pathways
# Recombinants of KOs from multiple individuals
# So, I can do an individual-by-individual run
# and isolate the 6 pathways, creating counts of individuals
# So I can have more confidence in the outputs
# I also toyed with the idea of running the backgound as the foreground...


# libraries ---------------------------------------------------------------
#libraries
library(duckdb)
library(tidyverse)
library(MicrobiomeProfiler)
library(ggplot2)


# bring in source object --------------------------------------------------------
# this will differ from Aarons method by using my custom background (pathway_analyser_comparatives.R)
# so the patterns in that and this could be wildly different
# but I stand by my method being right for adding relevance

# 1. load full object from local database
# i. open connection
con <- dbConnect(duckdb(), here::here("02_analysis/analysis_database/pseudomonas.db"), read_only=TRUE)

# ii. bring in object
pseudomonas_annotations <- dbGetQuery(con, 
                                      "SELECT * 
                          FROM pseudomonas_annotations pa;")

# iii. get custom background universe
unique_universe_kos <- dbGetQuery(con, 
                                  "SELECT * 
                          FROM pseudomonas_annotations pa;")
# same object, but I am just going to role with it for simplicity
# as this keeps the fore and backgrounds in seperate objects

# iv. close connection to duckdb
dbDisconnect(con)

# 1.2. get unique list of accessions to use in loop
accession_list <- unique(pseudomonas_annotations$Accession)


# 2. run enrichKO in a for loop so it enriches per individual -----------

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
  value <- enrichKO(gene = pseudomonas_annotations[pseudomonas_annotations$Accession == accession, "KEGG_ko"], universe = unique_universe_kos$KEGG_ko)
  # add the value to the list
  enriched_individuals <- append(enriched_individuals, values = value)
  # modify the name so I know what value relates to what accession
  names(enriched_individuals)[i] <- name
  i <- i + 1
}

# write object out to RDS so I dont have to rerun
saveRDS(enriched_individuals, file = "02_analysis/post_annotation_analysis/enrich_by_individual_myway_retry.rds")


# 3. make output into a dataframe -----------------------------------------
# easier to work with, good idea

enriched_individuals <- readRDS("02_analysis/post_annotation_analysis/enrich_by_individual_myway_retry.rds")
### create a dataframe of the key information ###

# 2 parts, the host info and the enrich info, join on accession?
# host info
# pseudo_info <- unique(data.frame("Host_group" = pseudomonas_annotations$animal_group,
#                                  "Accession" = pseudomonas_annotations$Accession))
# why do I create this? I dont use it, prob can be cut


# enrich info
y <- data.frame()

for (individual in names(enriched_individuals)) {
  z <- enriched_individuals[[individual]]@result
  z$Accession <- individual
  y <- rbind(y, z)
}

# why do I save it to RDS AND upload to duckdb?, I think I did some tinkering on duckdb
# so that can prob be ignored too
# write it out
saveRDS(y, file = "02_analysis/post_annotation_analysis/enriched_by_accession_myway_retry.rds")

# read out to duckdb
# con <- dbConnect(duckdb(), 
#                  dbdir = here::here("02_analysis/analysis_database/pseudomonas.db"),
#                  read_only = F)
# 
# duckdb::dbWriteTable(con, "enriched_by_accession", y)
# 
# dbDisconnect(con)

#y <- readRDS("02_analysis/post_annotation_analysis/enriched_by_accession_myway_retry.rds")

# 4. Tabulate / matrix results -----------------------------------------------------

# filter to just 6 pathways of interest
sig_pathway_comp <- y %>% 
  filter(ID %in% c("map00270", "map02040", "map01240", "map02025", "map03010", "map03070"))
# should be correct number, 6 pathways x 335 genomes = 2010

# add in group data for grouping
pseudo_info <- unique(data.frame("Host_group" = pseudomonas_annotations$animal_group,
                                  "Accession" = pseudomonas_annotations$Accession))
# add group data on
sig_pathway_comp <- inner_join(sig_pathway_comp, pseudo_info, by = join_by(Accession == Accession))

# NOTE: dont want to filter out non-significant pathways, as need them for total to make proportion

# TEST: new column of significance
# an attempt to keep it in 1 object 
sig_pathway_comp$significant <- ifelse(sig_pathway_comp$p.adjust < 0.05, "Significant", "Insignificant")
# NOTE: I cant believe I did that in only 2 tries

# group and count
sig_pathway_count <- sig_pathway_comp %>% 
  group_by(Host_group) %>% 
  count(ID, significant, name = "Count") %>% 
  ungroup() %>% 
  filter(significant == "Significant")

# total
sig_pathway_total <- sig_pathway_comp %>% 
  group_by(Host_group) %>% 
  count(ID, name = "Total") %>% 
  ungroup()

# join together 
sig_pathways <- left_join(sig_pathway_count, sig_pathway_total, by = c("Host_group" = "Host_group", "ID" = "ID"))
# 1st T'ed this one, with some stack overflow help

# make a proportion column
sig_pathways$proportion <- sig_pathways$Count / sig_pathways$Total


# make into matrix (very much like a heatmap) -----------------------------
# adapt heatmap code

# add colour so text distinguishable from background
sig_pathways$text_colour = ifelse(sig_pathways$proportion < 0.4, "white", "black")
# found a better way of doing it
# heat_source_long <- heat_source_long %>%
#   mutate(text_color = ifelse(proportion < 0.4, "white", "black"))

# basic plot
base <- ggplot(sig_pathways, aes(ID, Host_group, fill= proportion)) + 
  geom_tile() +
  scale_fill_viridis_c(direction = 1) +
  theme_minimal()

# add border to significant cells (where p < 0.05)
# got this off of google ai, realise I probably shouldnt have, but it works
developed <- base +
  # Add border to cells > threshold
  geom_tile(data = filter(sig_pathways, proportion > 0.8), 
            color = "black", size = 1, fill = NA) +
  labs(x = "KEGG Pathway ID", y = "Host Taxon", fill = "Proportion of\nindividuals") +
  geom_text(data = sig_pathways,
            aes(label = round(proportion, digits = 3), color = text_colour), size = 3) +
  scale_color_identity()
developed
# save
ggsave(filename = "03_outputs/02_enrichment/my_method_check_reran_again.jpeg", plot = developed, width = 17, height = 13, units = "cm")

# NOTE: Can upgrade with multiple confidence boxes, at like 80, 50 and 30. Like Claude.ai said would be good, rather than on/off
# may want to see whole output of this, instead of filtering to just the relevant pathways
# results are... funky bananas