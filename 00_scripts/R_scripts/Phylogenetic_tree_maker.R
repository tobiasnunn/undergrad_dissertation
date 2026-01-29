
# libraries ---------------------------------------------------------------

library(treeio)
library(here)
library(ape)
library(duckdb)
library(ggtree)
library(tidyverse)

# read in the whole tree and prune ----------------------------------------
# just want to keep tips in the list of accessions in the 4 groups

tree <- read.newick(here("02_analysis/gtdbtk_rerun/infer/gtdbtk.bac120.decorated.tree"))

#check the tree
print(tree)

# Get the tip labels
tip_labels <- tree$tip.label


# prune 1 - remove references ---------------------------------------------

# get list of my genomes
my_genome_tips <- tip_labels[grepl("_genome", tip_labels)]
length(my_genome_tips)

#prune tree
pruned_tree <- keep.tip(tree, my_genome_tips)

# Check result
print(paste("Pruned tree has", length(pruned_tree$tip.label), "tips"))

# now that the references are gone, the "_genome" is a hinderance to the second prune (and more generally the labelling)
# Remove "_genome" from the tip labels
pruned_tree$tip.label <- gsub("_genome$", "", pruned_tree$tip.label)

# check
head(pruned_tree$tip.label)


# prune 2 - remove samples not in 4 host groups ---------------------------

# bring in list of accessions in 4 host groups of interest from duckDB
# connect to duckdb 
con <- dbConnect(duckdb(), here::here("02_analysis/analysis_database/pseudomonas.db"), read_only=TRUE)

# get list of accessions with groups
study_pseudomonas <- dbGetQuery(con, 
                                  "SELECT DISTINCT pa.Accession as label, 
                                  pa.animal_group,
                                  pm.assemblyInfo.biosample.host,
                                  pm.organism.organismName, 
                                  pm.organism.infraspecificNames.strain,
                                  pm.assemblyInfo.biosample.isolationSource,
                                  pm.assemblyInfo.biosample.geoLocName,
                                  pm.assemblyInfo.biosample.latLon,
                                  pm.assemblyInfo.biosample.collectionDate
                                  FROM pseudomonas_metadata pm
                                  INNER JOIN pseudomonas_annotations pa 
                                  ON pa.Accession = pm.accession
                                  ORDER BY pm.organism.organismName;")
dbDisconnect(con)

# see how many species of Pseudomonas there are
unique(study_pseudomonas$organismName)
# Serpens galinarum is in the genus Pseudomonas apparently


# too many species to colour all, find top few species
study_pseudomonas %>% select(organismName) %>% group_by(organismName) %>% 
  count() %>% ungroup() %>% arrange(desc(n))


#match Pseudomonas species to colours
spec_colors <- c("Pseudomonas aeruginosa" = "#742881",
                 "Pseudomonas plecoglossicida" = "#986EAC",
                 "Pseudomonas putida" = "#E5D4E8",
                 "Pseudomonas chlororaphis" = "#D9F1D5",
                 "Pseudomonas fluorescens" = "#5CAE63",
                 "Pseudomonas guariconensis" = "#1B7939"
)

# add to study pseudomonas
study_pseudomonas <- study_pseudomonas %>% 
  mutate(big_species = ifelse(organismName %in% names(spec_colors), organismName, NA))


# keep only where tip.label is in study_pseudomonas
pruned_tree2 <- keep.tip(pruned_tree, study_pseudomonas$label)



# formatting and result ---------------------------------------------------
#match groups to colours
group_colors <- c(
  "Amphibian" = "#F0E442",
  "Bird" = "#009E73",
  "Fish" = "#56B4E9",
  "Mammal" = "#CC79A7")



## try both with pruned long branches --------------------------------------

#prune long branches
pruned_tree2$edge.length
sort(pruned_tree2$edge.length)
# drop last tip (which should match to the longest edge)
pruned_tree3 <- drop.tip(pruned_tree2, "GCA_052920985.1")


# tree -------------------------------------------------

# test plot (most basic)
tree_plot <- ggtree(pruned_tree3, layout = "circular") #roundrect or circular?
tree_plot 


# add on host information
tree_plot <- tree_plot %<+% study_pseudomonas

# building on what is above
colour_tree_plot <- tree_plot + 
  geom_tiplab(align = TRUE,
              linesize = 0.5,
              geom = "label", 
              aes(fill = animal_group, colour = animal_group),
              size = 1,
              label.padding = unit(0.35, "lines")) +
  scale_color_manual(name = "Animal Host Group", values = group_colors, guide = "none") +
  scale_fill_manual(name = "Animal Host Group", values = group_colors) +
  ggnewscale::new_scale_colour() +
  geom_segment(data = . %>% filter(isTip),
               aes(x = 0.35, xend = 0.375, y = y, yend = y, color = big_species), 
               size = 3) +
  scale_colour_manual(name = "Species", values = spec_colors, na.value = "lightgrey") +
  theme(legend.position = "right",
        legend.box = "vertical",  
        legend.spacing.y = unit(0.5, "cm"))

colour_tree_plot
#NOTES FOR TOMORROW, THINGS TO TEST / FIX:
# 1. fix the colour issue of the 1st ring DONE
# 2. see if viridis could mean I could get a colour for every species
# 3. maybe a third ring on another data attribute, like date or something
# 4. explore "ggtree extra"
# maybe try make it so strains are added to their species, could also just say "in many cases, the strain was included in the organismName field, so those did not make it to the species ring

# output to image file
ggsave(filename = "03_outputs/01_phylo/aligned_tree_2ring.jpeg", plot = colour_tree_plot, width = 27, height = 22, units = "cm")




# test area for latest tree

colour_tree_plot <- tree_plot + 
  geom_tiplab(align = TRUE,
              linesize = 0.5, 
              aes(colour = animal_group),
              size = 0.001) +
  geom_segment(data = . %>% filter(isTip),
               aes(x = 0.290, xend = 0.315, y = y, yend = y, color = animal_group), 
               size = 2) +
  scale_color_manual(name = "Animal Host Group", values = group_colors) +
  ggnewscale::new_scale_colour() +
  geom_segment(data = . %>% filter(isTip),
               aes(x = 0.320, xend = 0.345, y = y, yend = y, color = big_species), 
               size = 3) +
  scale_colour_manual(name = "Species", values = spec_colors, na.value = "lightgrey") +
  theme(legend.position = "right",
        legend.box = "vertical",  
        legend.spacing.y = unit(0.5, "cm"))

colour_tree_plot

# output to image file
ggsave(filename = "03_outputs/01_phylo/aligned_tree_2ring2.jpeg", plot = colour_tree_plot, width = 27, height = 22, units = "cm")


# add country info --------------------------------------------------------
# look for biases based on sampling site

# how many unique values
unique(study_pseudomonas$Country)


# too many species to colour all, find top few countries
study_pseudomonas %>% select(Country) %>% group_by(Country) %>% 
  count() %>% ungroup() %>% arrange(desc(n))

# mutate onto study_pseudomonas
study_pseudomonas <- study_pseudomonas %>% 
  separate_wider_delim(geoLocName, delim = ":", too_few = "align_start", names = c("Country", "region"))

# add continent info --------------------------------------------------------
# too many countries for ring to make sense, so going to continent-scale, other than China as it has 100+ samples in itself (i.e. I only really want to know if china is biasing the data)
# got Claude to sort them into groups, like I did with hosts

# Create the country-to-region mapping
country_mapping <- tibble::tribble(
  ~country,              ~wider_region,
  "China",               "China",  # Kept separate due to sample size (n=103)
  "Japan",               "Other Asia",
  "South Korea",         "Other Asia",
  "India",               "Other Asia",
  "Bangladesh",          "Other Asia",
  "Saudi Arabia",        "Other Asia",
  "Malaysia",            "Other Asia",
  "Turkey",              "Other Asia",  # Transcontinental, but primarily Asian
  
  "USA",                 "North America",
  "Canada",              "North America",
  "Mexico",              "North America",
  
  "United Kingdom",      "Europe",
  "France",              "Europe",
  "Czech Republic",      "Europe",
  "Portugal",            "Europe",
  "Spain",               "Europe",
  "Germany",             "Europe",
  "Switzerland",         "Europe",
  "Estonia",             "Europe",
  
  "Brazil",              "South America",
  
  "Australia",           "Oceania",
  
  "Tunisia",             "Africa",
  
  "Antarctica",          "Antarctica",
  "North Sea",           "Marine (unassigned)",  # Not a country
  
  NA_character_,         "Unknown",
  "not collected",       "Unknown"
)

# Convert to factor with explicit level ordering for your tree
region_order <- c("China", "Other Asia", "North America", "Europe", 
                  "South America", "Oceania", "Africa", "Antarctica", 
                  "Marine (unassigned)", "Unknown")

country_mapping <- country_mapping %>%
  mutate(wider_region = factor(wider_region, levels = region_order))

# add that info onto study pseudomonas
study_pseudomonas <- study_pseudomonas %>% 
  inner_join(country_mapping, by = join_by(Country == y$country))


# rerun tree_plot with new study_pseudomonas
tree_plot <- ggtree(pruned_tree3, layout = "circular") #roundrect or circular?
tree_plot 


# add on host information
tree_plot <- tree_plot %<+% study_pseudomonas


# add third ring for continent (ideally using scale colour viridis)
colour_tree_plot <- tree_plot + 
  geom_tiplab(align = TRUE,
              linesize = 0.5, 
              aes(colour = animal_group),
              size = 0.001) +
  geom_segment(data = . %>% filter(isTip),
               aes(x = 0.290, xend = 0.315, y = y, yend = y, color = animal_group), 
               size = 2) +
  scale_color_manual(name = "Animal Host Group (inner ring)", values = group_colors,
                     guide = guide_legend(order = 3)) +
  ggnewscale::new_scale_colour() +
  geom_segment(data = . %>% filter(isTip),
               aes(x = 0.320, xend = 0.345, y = y, yend = y, color = big_species), 
               size = 3) +
  scale_colour_manual(name = "Species (middle ring)", 
                      values = spec_colors, na.value = "lightgrey",
                      guide = guide_legend(order = 2)) +
  ggnewscale::new_scale_colour() +
  geom_segment(data = . %>% filter(isTip),
               aes(x = 0.350, xend = 0.375, y = y, yend = y, color = wider_region), 
               size = 3) +
  scale_color_viridis_d(name = "Collection Region (outer ring)", direction = -1,
                        guide = guide_legend(order = 1)) +
  theme(legend.position = "right",
        legend.box = "vertical",  
        legend.spacing.y = unit(0.5, "cm"))

colour_tree_plot
# output to image file
ggsave(filename = "03_outputs/01_phylo/aligned_tree_3ring.jpeg", plot = colour_tree_plot, width = 27, height = 22, units = "cm")


# update middle ring OPTIONAL ------------------------------------------------------
# The main point of showing that p. aeruginosa is huge and links with mammals is shown
# but maybe if I remove strain from some of them some of the groups will increase in size a bit

# mutate onto study_pseudomonas
study_pseudomonas_copy <- study_pseudomonas %>% 
  separate_wider_delim(organismName, delim = " ", too_many = "debug", names = c("Genus", "species")) %>% select(-c(organismName_ok, organismName_pieces, organismName))

study_pseudomonas_copy$binomial <- paste(study_pseudomonas_copy$Genus,
                                         study_pseudomonas_copy$species)
  
study_pseudomonas_copy <- study_pseudomonas_copy %>% select(c(
  label, animal_group, host, binomial, organismName_remainder, strain, big_species, 
  wider_region, Country, region, latLon, isolationSource, collectionDate))

# rerun the stuff that makes "big_species"
# see how many species of Pseudomonas there are
unique(study_pseudomonas_copy$binomial)
# Serpens galinarum is in the genus Pseudomonas apparently


# too many species to colour all, find top few species
study_pseudomonas_copy %>% select(binomial) %>% group_by(binomial) %>% 
  count() %>% ungroup() %>% arrange(desc(n))

#see if Pseudomonas sp. with the strain can be used
sp_tester <- filter(study_pseudomonas_copy, binomial == "Pseudomonas sp.") %>% 
  select(binomial, strain)
nrow(unique(sp_tester))
# right 76, so no, sp. with strain does not yield anything

# CONCLUSION: NOT WORTH UPDATING TREE FOR THIS INFORMATION

# statistical testing for bias / confounding ------------------------------

# as I thought, the chi squared
# Claude also suggested Cramér's V

# (for now) DECIDED NOT TO DO AS DESCRIPTIVE STATS ARE LIKELY MORE INFORMATIVE IN THIS CASE (i.e. the tree shows 140/ 170 mammal samples belong to p. aeruginosa (~80%), showing taxonomic bias.)


########################### DEVO stuff ##################################
## circle (no length) ------------------------------------------------------

# cross-examine in the light of phylogenetic information ------------------

# bring in / create an object containing bacterial taxonomic information

### Things to think of:-----------------
# highlighting - done here
# gradients - would be fun
# tiplab2 - for circular stuff
# tiny segment - not sure what I meant by that
# shiny - interactive

# test plot (most basic)
tree_plot_circ <- ggtree(pruned_tree2, layout = "circular", branch.length = "none")
tree_plot_circ 




# add on host information
tree_plot_circ <- tree_plot_circ %<+% study_pseudomonas

# add colour
colour_tree_plot <- tree_plot_circ + 
  geom_tippoint(aes(color = animal_group), size = 1) +
  geom_segment(
    data = . %>% filter(isTip),  # Only apply to tips, not internal nodes
    aes(x = x + 3, xend = x + 15, y = y, yend = y, color = animal_group),
    linewidth = 2.5, alpha = 1
  ) +
  scale_color_manual(values = group_colors) +  
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold")) +
  labs(color = "Animal Host Group")
colour_tree_plot


# rect (inc length) -------------------------------------------------------

# test plot (most basic)
tree_plot_rect <- ggtree(pruned_tree2, layout = "roundrect")
tree_plot_rect 




# add on host information
tree_plot_rect <- tree_plot_rect %<+% study_pseudomonas

# add colour
colour_tree_plot <- tree_plot_rect + 
  geom_tippoint(aes(color = animal_group), size = 1) +
  geom_segment(
    data = . %>% filter(isTip),  # Only apply to tips, not internal nodes
    aes(x = x + (0.1 * x), xend = x + (0.33 * x), y = y, yend = y, color = animal_group),
    linewidth = 1.5, alpha = 0.7
  ) +
  scale_color_manual(values = group_colors) +  
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold")) +
  labs(color = "Animal Host Group")
colour_tree_plot

# test plot (most basic)
tree_plot <- ggtree(pruned_tree3, layout = "roundrect") #roundrect or circle?
tree_plot 


# add on host information
tree_plot <- tree_plot %<+% study_pseudomonas

# add colour
colour_tree_plot <- tree_plot + 
  geom_tippoint(aes(color = animal_group), size = 1) +
  geom_segment(
    data = . %>% filter(isTip),  # Only apply to tips, not internal nodes
    aes(x = x + 0.001, xend = x - 0.05, y = y, yend = y, color = animal_group),
    linewidth = 1, alpha = 0.7
  ) +
  scale_color_manual(values = group_colors) +  
  theme(legend.position = "right",
        legend.title = element_text(face = "bold")) +
  labs(color = "Animal Host Group") 
colour_tree_plot
# highlighting
# gradients
# tiplab2