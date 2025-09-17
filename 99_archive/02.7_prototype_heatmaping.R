
# libraries ---------------------------------------------------------------

library(tidyverse)  # all purpose tool
library(RPostgres)  # allows issueing of commands to postgres database
library(DBI)        # controls connection to database
library(ggplot2)    # good graphics
library(scales)     # allows the blocking for styling and organisation
 
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

KO_files <- data.frame(accession = list.files(path = "../01_inputs/02_previous_analysis_results/eggnog-mapper/", pattern = "*ko_pathway.tsv"),
                       file_path = list.files(path = "../01_inputs/02_previous_analysis_results/eggnog-mapper/", pattern = "*ko_pathway.tsv", full.names = T))
KO_files$accession <- gsub("_ko_pathway.tsv", "", KO_files$accession)
KO_files$accession_number <- gsub("GCA_|GCF_|\\.[0-9]+$", "", KO_files$accession) # modify the gsub so that any number after the point at the end gets cut


# as part of bringing the file in I have to manually add the local samples because
# the file was created purely from NCBI data
accession_with_host <- read_delim("../03_outputs/host_evaluation/accession_to_host_class.tsv") %>% 
  mutate(accession_number = gsub("GCA_|GCF_|\\.[0-9]+$", "", accession)) %>% 
  add_row(host_value = "dendrobates tinctorius", class_name = "Amphibia", class_tax_id = 8292, genus_id = 33882, genus_name = "Microbacterium",
          accession = "2Dt1l", accession_number = "2Dt1l") %>% 
  add_row(host_value = "dendrobates tinctorius", class_name = "Amphibia", class_tax_id = 8292, genus_id = 33882, genus_name = "Microbacterium",
          accession = "3Dt2h", accession_number = "3Dt2h") %>% 
  filter(accession_number %in% KO_files$accession_number)
#TODO: at some point I will need to do the other 8 (im going to need to figure out that 1 sample I have yet to)

#get names of top 5 classes (these are the ones Microbacterium has right now)
top_6_classes <- accession_with_host %>% 
  group_by(class_name, genus_name) %>% 
  count(name = "count_of_accessions") %>% 
  pivot_wider(names_from = genus_name, values_from = count_of_accessions, names_sort = TRUE) %>% 
  rowwise() %>% 
  mutate(Total = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  ungroup()  %>%
  arrange(desc(Total)) %>% 
  filter(Microbacterium>1)

# filter 
filtered_accessions <- accession_with_host %>% 
  filter(genus_name == "Microbacterium") %>% 
  filter(class_name %in% top_6_classes$class_name) %>% 
  left_join(KO_files, by=join_by(accession_number==accession_number))



# aggregation
mapids_micro <- data.frame(accession = character(), class = character(), map_id = character(), map_description = character(),
                     gene_ratio = character(), p_val = numeric())

for (i in 1:nrow(filtered_accessions)) {
  #i <- 1
  filename <- filtered_accessions$file_path[i]
  if(file.exists(filename)){
    intermed <- read_delim(filename, col_names = TRUE, delim = "\t")
    accession <- filtered_accessions$accession.x[i]
    class <- filtered_accessions$class_name[i]
    map_id <- intermed$ID
    map_description <- intermed$Description
    gene_ratio <- intermed$GeneRatio
    p_val <- intermed$pvalue
    mapids_micro <- add_row(mapids_micro, accession = accession, class = class, map_id = map_id, map_description = map_description,
                      gene_ratio = gene_ratio, p_val = p_val)
  }
}


# outputs
mapids_summary_absolute_byclass <- mapids_micro %>% 
  group_by(class) %>% 
  count(map_id, name = "count") %>% 
  pivot_wider(names_from = class, values_from = count, values_fill = 0)

# proportionising
# lowered hi lim and raised lo lim to get more out of small dataset
high_lim <- 0.6
low_lim <- 0.5

# cells need to be proportioned based on group total, this frame gets those totals
class_count <- filtered_accessions %>% group_by(class_name) %>% count(name = "total")

#TODO: at some point, I will need to modify the code to do cutdown heatmaps
# on the event that one group is over-represented, however, the groups are on
# a much more even evolutionary distance than when I was doing genera.

#genera <- c("Brevibacterium", "Microbacterium", "Brachybacterium", "Pantoea", "Sphingomonas")

#raw_counts <- select(raw_counts, map_id | all_of(genera))

proportional_mapids <- mapids_summary_absolute_byclass %>% 
  pivot_longer(cols = -map_id, names_to = "class_name", values_to = "count") %>% 
  left_join(class_count, by = "class_name") %>% 
  mutate(prop = count / total) %>% 
  select(map_id, class_name, prop) %>% 
  pivot_wider(names_from = "class_name", values_from = "prop") 

numeric_cols <- names(proportional_mapids)[sapply(proportional_mapids, is.numeric)]

filtered_prop_count <- proportional_mapids %>% 
  rowwise() %>%
  filter(sum(c_across(where(is.numeric)) > high_lim) == 1 & 
           sum(c_across(where(is.numeric)) < low_lim) == length(numeric_cols) - 1) %>%
  ungroup() %>%
  pivot_longer(cols = -map_id, names_to = "class", values_to = "prop")

# metadata for boxes
kegg_metadata <- read_delim("../01_inputs/02_previous_analysis_results/eggnog-mapper/genera_KEGG_metadata.tsv", delim = "\t") %>% 
  rename(kegg_class = class)

filtered_prop_count <- filtered_prop_count %>% 
  left_join(kegg_metadata, by = join_by(map_id == map)) %>% 
  select(map_id, class, prop, name, kegg_class, subclass) %>% 
  mutate(yellow = prop > 0.5)

filtered_prop_count$perc <- scales::percent(filtered_prop_count$prop, accuracy = 0.1)

# heatmap
proto_class_heatmap <- ggplot(data = filtered_prop_count, mapping = aes(x = fct_rev(class),
                                                                   y = name, 
                                                                   fill = prop)) +
  geom_tile(colour = "lightgrey", lwd = 0.5, linetype = 1) +
  labs(x =  "Host Class", y ="KEGG Pathway", fill = "proportion\nenriched\ngenomes",
       title = "PROTOTYPE\nComparative prevalance of KEGG pathways between\ntaxonomic Class of host for bacterial genus\nMicrobacterium",
       subtitle = "n = 49 (10, 28, 2, 4, 3, 2), limits for inclusion 60 and 50%") + # there has to be a way to automate that
  theme(axis.text.x = element_text(vjust = 0.5), text = element_text(size = 14)) +
  #facet_grid(subclass ~ ., scales = "free", space = "free") +
  scale_fill_viridis_c(limits = c(0,1), option = "plasma") +
  theme(strip.placement = "outside") +
  theme(strip.text.y = element_text(angle = 0), strip.text = element_text(size = 16)) +
  theme(axis.text = element_text(size = 10), plot.title = element_text(size = 16)) +
  facet_wrap(~subclass, scales="free_y", nrow = 5, ncol = 1) + 
  geom_text(aes(label = perc, colour = yellow)) +
  scale_colour_manual(values=c("white", "black")) + 
  guides(colour = "none")

proto_class_heatmap
# scales="free_y" is instrumental in getting the x-axis to not repeat

ggsave("../03_outputs/prototypes/proto_class_heatmap.png", 
       plot = proto_class_heatmap, width = 3000, height = 2500, units = "px")
