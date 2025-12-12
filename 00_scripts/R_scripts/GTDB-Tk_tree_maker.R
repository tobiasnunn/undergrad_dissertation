library(ggtree)
library(treeio)
library(tidyverse)

# Note: to save on time, in the file "gtdbtk.bac120.decorated.tree" the outgroup 
# was cut using the tool dendroscope
# following code from book by G Yu: https://yulab-smu.top/treedata-book/chapter2.html
# read in outputs of gtdbtk -----------------------------------------------
#shapes <- c(1:15)

tree <- treeio::read.newick(here::here("02_analysis/02_gtdbtk/pruned_gtdbtk.bac120.decorated.tree"))
tree_plot <- ggtree(tree, layout = "circular", branch.length = "none") %<+% mapping_data + 
  aes(color=host_group) +
  geom_tippoint(aes(shape = host_group, color = host_group))
tree_plot 
ggsave(here::here("02_analysis/02_gtdbtk/innitial_tree.png"), tree_plot, 
                  width = 20, height = 15)
#+ geom_text(aes(x=branch, label=tree$tip.label, vjust=-.5))
#?ggtree
#print(tree %>% as.treedata %>% as_tibble, n=30)
#print(bazinga@data)

# read in the metadata
groups <- readRDS(here::here("01_inputs/03_metadata/filtered_accessionframe.rds"))


mapping_data <- filter(tree_plot@data, tree_plot@data$isTip == TRUE) %>% 
  mutate(join_column = gsub("_genome", "", label)) %>% 
  left_join(groups, by = join_by(join_column == accession_id)) %>%
  mutate(
    megagroup = if_else(
      host_group %in% c("bovine_group", "canine_group"),
      host_group,
      NA_character_
    )
  )

p <- ggtree(tree, layout = "circular", branch.length = "none") %<+% mapping_data + 
  aes(color=megagroup)
p
