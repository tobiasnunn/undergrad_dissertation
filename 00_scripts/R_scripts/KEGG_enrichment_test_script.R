if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

BiocManager::install("GO.db")

library(clusterProfiler)
library(tidyverse)

# test code from book
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
kj <- kk@result
head(kk)

KO_file <- read_tsv("C:/Users/tobyn/OneDrive/Work/tnunn_research/2025 file storage/annotation_200_accessions.emapper.annotations.copy.gz")
KO_file$accession <- gsub("_[0-9]+$", "", KO_file$`#query`)
first_sample <- filter(KO_file, accession == "GCA_016429035.1")
unique(KO_file$accession)

# get just list of KOs
first_KO_list <- data.frame(KOs = first_sample$KEGG_ko) %>% 
  mutate(KOs = gsub("ko:", "", KOs)) %>% 
  separate_longer_delim(KOs, ",") %>% 
  filter(!KOs == "-")

# test code from book with my KOs
?enrichKEGG
first_result <- enrichKEGG(first_KO_list$KOs,
                           organism     = 'ko',
                           pvalueCutoff = 0.05)
first_resres <- first_result@result
head(first_resres)


# test inbuilt visualisations
library(enrichplot)  # Load visualization extensions

# 1. Bar plot - shows top enriched pathways
barplot(first_result, showCategory = 20)
?barplot.enrichResult
# just gives a basic plot of size against pathway, nothing crazy

# 2. Dot plot - combines significance and gene ratio
dotplot(first_result, showCategory = 20)
# bit more useful, showing the GeneRatio, as well as the absolute size

# 3. Gene-Concept Network - shows which genes belong to which pathways
cnetplot(first_result, categorySize = "pvalue", foldChange = NULL)
?cnetplot.enrichResult
# definitely somehting interesting, but maybe dont do on whole dataset
# it shows which KOs are in what pathways, making a fun-funky map, def
# something to put in the dissertation for wow factor

# 4. Heatmap plot - if you have multiple samples
# heatplot(enrichkegg_result)
# havent KEGG-ed multiple samples yet, would be good to do in future

# 5. Enrichment Map - shows pathway relationships
emapplot(first_result)
?`emapplot,enrichResult-method`
# not sure how to make it work, think it is very old and not still maintained

# 6. Ridge plot - shows gene distribution across pathways
ridgeplot(first_result)
?ridgeplot.gseaResult
# cant do on enrich result

# 7. Tree plot - hierarchical view of pathways
treeplot(first_result)
?treeplot.enrichResult
?pairwise_termsim.enrichResult
# "similarity matrix?" - apparently needed to do 5 and 7
first_result_pairwise <- pairwise_termsim(first_result)
treeplot(first_result_pairwise)
# no idea why this now works, but it does
# produces an interesting, if freaky diagram

# 5. Enrichment Map - shows pathway relationships
emapplot(first_result)
?`emapplot,enrichResult-method`
# not sure how to make it work, think it is very old and not still maintained
emapplot(first_result_pairwise)
# also needs the freaky thing to work
# looks to be like 3, but at the pathway level, so simpler, but less info