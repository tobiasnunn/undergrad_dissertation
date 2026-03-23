# libraries
library(flextable)
library(duckdb)
library(tidyverse)

# table 1 - distribution of bacteria between host groups----------------------------------------------

con <- dbConnect(duckdb(), here::here("02_analysis/analysis_database/pseudomonas.db"), read_only=TRUE)

# table of host counts
host_group_counts <- dbGetQuery(con, 
                                "SELECT animal_group as 'Animal group', COUNT(*) as Count
FROM (SELECT DISTINCT Accession, animal_group
FROM pseudomonas_annotations pa) as Q1
GROUP BY animal_group;")

# close connection to duckdb
dbDisconnect(con)

# order biggest to smallest
host_group <- host_group_counts[order(-host_group_counts$Count),]
# host_group_dplyr <- host_group_counts %>% arrange(desc(Count))

# scientific names as they will want that
host_group$`Animal group` <- c("Mammal (Mammalia)", "Fish (Osteichthyes)",
                                 "Bird (Aves)", "Amphibian (Amphibia)")

# flextable
table1 <- flextable(host_group) %>% 
  set_header_labels(`Animal group` = "Host group") #%>%
  #border_inner_v(part="all", border = fp_border_default())

table1

save_as_image(table1, "03_outputs/03_descriptive_tables/table1.png")  


# table 2 - breakdown by species ------------------------------------------

# right, I will have done something like this when making the tree, maybe I can
# pull code from there

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
table2 <- study_pseudomonas %>% select(organismName) %>% group_by(organismName) %>% 
  count() %>% ungroup() %>% arrange(desc(n))
# check to make sure that sums to 337
#sum(table2$n)
#yep, we good

# too many rows, remember to note there are more
table2 <- head(table2, 10)

# flextable
flextable2 <- flextable(table2)

flextable2
save_as_image(flextable2, "03_outputs/03_descriptive_tables/table2.png")  


# table 3 - KO lists of groups vs background ------------------------------

# hmm, how am I going to make this one?
# pathway in group | universe pathways | percent ?

con <- dbConnect(duckdb(), here::here("02_analysis/analysis_database/pseudomonas.db"), read_only=TRUE)

# bring in count and percentage of unique lists of KOs
unique_KO_proportions <- dbGetQuery(con, 
                                  "WITH total_distinct_ko AS 
(SELECT COUNT(*) AS total_ko FROM (SELECT DISTINCT KEGG_ko
FROM pseudomonas.main.pseudomonas_annotations)
),
count_by_group AS 
(SELECT animal_group, COUNT(*) AS total_by_hosts FROM (SELECT DISTINCT animal_group, KEGG_ko
FROM pseudomonas.main.pseudomonas_annotations)
GROUP BY animal_group
UNION 
SELECT 'Universe' AS animal_group, COUNT(*) AS total_ko FROM (SELECT DISTINCT KEGG_ko
FROM pseudomonas.main.pseudomonas_annotations))
SELECT cbg.animal_group, cbg.total_by_hosts, tdk.total_ko, cbg.total_by_hosts/tdk.total_ko AS percent_of_total
FROM count_by_group cbg, total_distinct_ko tdk
ORDER BY percent_of_total DESC;")


# close connection to duckdb
dbDisconnect(con)

# flextable
flextable3 <- flextable(select(unique_KO_proportions, -total_ko)) %>% 
  set_header_labels(animal_group = "Host group", total_by_hosts = "Distinct KO by group",
                    percent_of_total = "Proportion of universe KOs") %>% 
  mk_par(j = 3, part = "body",
    value = as_paragraph(as_chunk(percent_of_total, formatter = fmt_pct))) %>% 
  width(width = 1.5) %>% 
  bg(i = 1, bg = "#ededed")

#save as img
save_as_image(flextable3, "03_outputs/03_descriptive_tables/table3.png")  
