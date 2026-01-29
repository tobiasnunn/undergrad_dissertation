
# libraries ---------------------------------------------------------------

library(ggplot2)
library(plo)

# load in output of pathway analysis --------------------------------------

diff_paths <- readRDS("02_analysis/post_annotation_analysis/differential_pathways.rds")



# heatmap -----------------------------------------------------------------

# create the base heatmap source item
heatmap_source <- diff_paths %>%  select(c(ID, p.adjust.amphibian,
                                           p.adjust.mammal, p.adjust.bird,
                                           p.adjust.fish))

# make it long
heat_source_long <- heatmap_source %>% 
  pivot_longer(cols = c(p.adjust.amphibian,
                        p.adjust.mammal, p.adjust.bird,
                        p.adjust.fish), names_to = "p.adjusted", 
               values_to = "p_value")

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
  labs(x = "KEGG Pathway ID", y = "adjusted P value")

# save
ggsave(filename = "03_outputs/01_phylo/diff_paths_heatmap.jpeg", plot = developed, width = 20, height = 15, units = "cm")




# extra stuff if there is time --------------------------------------------


# interactivity with plotly
# having an error: Error in unpackPkgZip(foundpkgs[okp, 2L], foundpkgs[okp, 1L], lib, libs_only,  : 
# ERROR: failed to lock directory ‘C:\Users\tobyn\AppData\Local\R\win-library\4.5’ for modifying
#Try removing ‘C:\Users\tobyn\AppData\Local\R\win-library\4.5/00LOCK’


# basically, add a "text" column onto the end of heat_source_long, with the stuff in it formatted in the way the text should be, then call that in the ggplotly()