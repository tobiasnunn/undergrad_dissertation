
# 01. libraries -----------------------------------------------------------


library(tidyverse)
library(ggplot2)
library(flextable)
library(scales)

# 02. read in metadata files ----------------------------------------------


genus_tax_info <- readRDS(here::here("01_inputs/03_metadata/taxonomy_information_Sphingomonas_Brevibacterium_Microbacterium_Brachybacterium_Pantoea_2025-09-10.rds"))
accession_info_by_genus <- readRDS(here::here("01_inputs/03_metadata/accession_information_43668_1696_33882_53335_13687_2025-09-10.rds"))


# 03. extract a list of host information ----------------------------------

## MK 1 - base R but imperfect because uses coalesce from purrr
# bazinga <- data.frame(genus = as.character(), accession = as.character(), host = as.character(), iso_source = as.character())
# 
# 
# for (i in seq_along(accession_info_by_genus)) {
#   for (j in seq_along(accession_info_by_genus[[i]])) {
#     bazinga <- rbind(bazinga, data.frame(genus = names(accession_info_by_genus[i]), 
#                                          accession = names(accession_info_by_genus[[i]][j]),
#                                          host = coalesce(accession_info_by_genus[[i]][[j]]$assembly_info$biosample$host, NA),
#                                          iso_source = coalesce(accession_info_by_genus[[i]][[j]]$assembly_info$biosample$isolation_source, NA)))
#   }
# }

## MK 2 - base R without additions (Claude assisted, I really should look at that page about R shortcuts)
# bazinga2 <- data.frame(genus = as.character(), accession = as.character(), host = as.character(), iso_source = as.character())
# 
# for (i in seq_along(accession_info_by_genus)) {
#   for (j in seq_along(accession_info_by_genus[[i]])) {
#     bazinga2 <- rbind(bazinga2, data.frame(genus = names(accession_info_by_genus[i]), 
#                                          accession = names(accession_info_by_genus[[i]][j]),
#                                          host = accession_info_by_genus[[i]][[j]]$assembly_info$biosample$host %||% NA,
#                                          iso_source = accession_info_by_genus[[i]][[j]]$assembly_info$biosample$isolation_source %||% NA))
#   }
# }

## MK 3 - tidyverse, but imperfect because I still have to use a lot of base R notation
# banana <- data.frame(genus = as.character(), accession = as.character(), host = as.character(), iso_source = as.character())
#   
# for (i in seq_along(accession_info_by_genus)) {
#   for (j in seq_along(accession_info_by_genus[[i]])) {
#     banana <- rbind(banana, data.frame(genus = names(pluck(accession_info_by_genus[i])),
#                                        accession = names(pluck(accession_info_by_genus[[i]][j])),
#                                        host = pluck(accession_info_by_genus, i, j, "assembly_info", "biosample", "host", .default = NA),
#                                        iso_source = pluck(accession_info_by_genus, i, j, "assembly_info", "biosample", "isolation_source", .default = NA)))
#     
#     }
# }

## MK 4 - tidyverse, without base R (Claude assisted, suppose I could read a paper about tidyverse...)
banana2 <- imap(accession_info_by_genus, function(taxon_data, tax_id){
  imap(taxon_data, function(accession_data, accession_id){
    data.frame(
      tax_id = as.integer(tax_id),
      accession_id = accession_id,
      host = pluck(accession_data, "assembly_info", "biosample", "host", .default = NA),
      isolation_source = pluck(accession_data, "assembly_info", "biosample", "isolation_source", .default = NA)
    )
  }) %>% list_rbind()
}) %>% list_rbind()


# 04. manipulate data for visualisation -----------------------------------


# bring in genus name
joined_data <- left_join(banana2, genus_tax_info, by="tax_id") %>% 
  select(organism_name, tax_id, accession_id, host, isolation_source)

# split into host and iso_source dataframes
host_data <- select(joined_data, -isolation_source) %>% 
  mutate(standardised_host = toupper(host))
iso_data <- select(joined_data, -host) %>% 
  mutate(standardised_iso = toupper(isolation_source))

# unique on standardised host name
unique_host_data <- data.frame(unique(host_data$standardised_host))
# save it to excel so I can parse it for stuff to 86
write.csv(unique_host_data, file = "02_analysis/metadata/unique_host_data.csv")
nix_values <- c('MISSING','MSSING','N/A','NA', NA, 'NO COLLECTED','NOT APPLICABLE','NOT AVAILABLE','NOT COLLECTED','NOT PROVIDED','UNKNOWN')

# create valid column to see what % has a host or iso_source
host_data$valid <- !host_data$standardised_host %in% nix_values


# 04.5 - further manipulation for more complex graphs ---------------------

valid_hosts <- filter(host_data, valid == TRUE) %>% 
  select(standardised_host) %>% 
  unique() %>% 
  arrange(standardised_host)

counted_hosts <- filter(host_data, valid == TRUE) %>% 
  select(-c(host,tax_id, valid)) %>% 
  count(organism_name, standardised_host) %>% 
  group_by(organism_name) %>% 
  slice_max(n , n = 5, with_ties = F) %>% 
  pivot_wider(names_from = organism_name, values_from = n)


# 05. visuals -------------------------------------------------------------


## 1. % with host
host_perc_frame <- host_data %>% count(valid) %>% 
  mutate(perc = n/sum(n)*100, label = paste0(n, "\n", round(perc, digits = 1), "%"))

host_perc <- ggplot(host_perc_frame, aes(x = valid, y = n)) +
  geom_col(fill = c("navy","lightblue")) +
  theme_minimal() +
  geom_text(aes(label = label, y = n/2, colour = ifelse(valid == FALSE, "white", "black")), size = 5) +
  scale_color_identity()
host_perc

## 2. breakdown of top 5 hosts per genus (matrix)

top5_per_genera <- flextable(counted_hosts) %>% 
  set_header_labels(standardised_host = "Host") %>% 
  align(i = NULL, j = ~. -standardised_host, align = "center", part = "body") %>% 
  color(color = "white", part = "body", j = ~. -standardised_host)  %>% 
  color(color = "black", part = "body", j = 5, i = 13) %>% 
  bg(j = ~. -standardised_host, bg = scales::col_numeric(palette = "plasma", 
                                                         domain = c(0, 145), na.color = "white"), part = "body")
top5_per_genera

# fancy Claude way of doing it so I dont have to hard code the font colour

# top5_per_genera2 <- flextable(counted_hosts) %>% 
#   set_header_labels(standardised_host = "Host") %>%
#   bg(j = ~. -standardised_host, bg = scales::col_numeric(palette = "plasma", 
#                                                          domain = c(0, 145), na.color = "white"), part = "body")
# 
# # Using column numbers
# for (j in 2:6) {
#   col_name <- names(counted_hosts)[j]
#   top5_per_genera2 <- top5_per_genera2 %>%
#     color(
#       i = as.formula(paste0("~ is.na(`", col_name, "`) | `", col_name, "` <= 80")),
#       j = j,
#       color = "white"
#     )
# }
# 
# top5_per_genera2
