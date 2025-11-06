file_name <- here::here("01_inputs/03_metadata/filtered_accessionframe.rds")
all_genera_accession_list <- readRDS(file_name)

#accession_list <- all_genera_accession_list[[as.character(taxon_ID)]]
accession_list <- all_genera_accession_list$accession_id
print(file_name)
print(nrow(all_genera_accession_list))
print(length(accession_list))

# make sure accessions are not being repeated
file_list <- data.frame(filename = list.files(path = here::here("01_inputs/04_fastas/"), pattern = "*.rds"))
file_list$accession <- gsub(".rds", "", file_list$filename)
accession_list <- accession_list[!accession_list %in% file_list$accession]
print(length(file_list$accession))
print(length(accession_list))

log_msg <- paste(Sys.time(), "Number of accessions to process:", length(accession_list))
#cat(log_msg, "\n", file = here::here(log_file), append = TRUE)
print(log_msg)
