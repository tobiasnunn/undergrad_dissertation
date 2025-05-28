
# libraries ---------------------------------------------------------------


library(readxl)     #reads in .xlsx formatted files
library(tidyverse)  # all purpose tool
library(httr2)      # used to access NCBI REST API
library(flextable)  # tablulation beautification
library(jsonlite)   # read in .json after API

# 01. collate autoMLST2 outputs -------------------------------------------


# as noted in notebook-week 2, 7 samples are in a .xlsx file, and 3 are in text files

main_data <- read_xlsx('../03_outputs/autoMLST2_subproject/other_7/Results from AutoMLST on Bangor samples.xlsx',
                       range = "B5:J105")  # reads in output table from screenscrape for 7/10 samples

other_files_names <- list.files(path = "../03_outputs/autoMLST2_subproject/", 
                                pattern = "^mash_distances.*\\.txt$", 
                                recursive = TRUE, full.names = TRUE)
other_files_names
# looks for files inside my github repo directory structure. It uses regex to 
#find the matching names (^ means the letter after it should be the first letter,
# the * means anything can be there, the $ means end of string). It is recursive,
# meaning it looks in all the subfolders in autoMLST2_subproject.

other_files <- lapply(other_files_names, read_delim) # reads named files into a list, one can note
# the matching names between the objects in the list and main_data. Also, these files are very large
# which helps explain why the other 7 didnt work, because theorhettically, that file was 7x larger

names(main_data) <- names(other_files[[1]]) # in order to join the objects, the column names
# need to be the same, this was not the case, so this line matches main_data to other_files

mlst_table <- bind_rows(other_files, main_data) # binds the rows of the two objects together
# the order of the two objects does not matter, "bind_rows" comes from dplyr in tidyverse


mlst_table$accession <- gsub(".*QS--flye_asm_(.*)_part2\\.fa.*", "\\1", mlst_table$`#Query_org`)
# In order to match with the checkm2 and gtdbtk data, I need a column of just the 
# accession name. This line uses very complex regex to extract the accession from
# inside the string in `#Query_org`.
                               
mlst_table <- mlst_table %>% group_by(accession) %>% slice_min(order_by = MASH_distance, n = 1) %>% 
  select(accession, Reference_assembly, Ref_name, MASH_distance, Estimated_ANI)
# takes mlst_table, groups by accession, I take only the best match by group.
# This is done with min MASH_distance, but could also be done with max Estimated_ANI
# because both are metrics of genetic closeness. I then keep only relevant columns

# I dropped columns "P-value", "Genus", "Order" and "Type_value" for space
# Type_value was interesting, it means if a reference assembly is a reference sample
# i.e. on the NCBI website they have a sticker saying "reference"

rm(main_data)
rm(other_files)
rm(other_files_names)
#cleaning

# 02. bring in gtdb-tk data -----------------------------------------------


gtdbtk <- read_delim('../01_inputs/previous_analysis_results/gtdb-tk/classify_wf/gtdbtk.bac120.summary.tsv')
# this .tsv is just an output of classify_wf, i did not do anything to it.
gtdbtk$accession <- gsub("flye_asm_(.*)_part2", "\\1", gtdbtk$user_genome)
# need to isolate the accession for combining

gtdbtk$closest_placement_taxonomy
gtdbtk$closest_placement_ani
gtdbtk$warnings
# not all of them have identifiers, I think because the MASH distance
# or some ANI value was set too high, so I might go back through and 
# change that value if I can / it exists to see what comes out


gtdbtk$closest_tax <- gsub(".*s__", "", gtdbtk$closest_placement_taxonomy)

gtdbtk <- gtdbtk %>%  select(
  accession, closest_tax, closest_placement_reference, closest_placement_ani,
  msa_percent
)
# need to cut a load of columns, otherwise it will bloat, so here are my comments
# on the stuff I have had to delete, for adding to notebook
# only 1Dt2d got a fastani value, due it being >95% ANI
# closest placement radius is 95.0
# "closest_placement_af" was cut because I dont know its relevance
    # had a look, it means "alignment fraction"(?)
# "classification_method" was cut, 9 were "topology and ANI", 1Dt100h was done using "RED"
# "note" was cut because it did not contain new info
# classify also provided more closely related samples, but that is not relevant to this
# translation table is 11 / other column lost was "red_value", no idea
# "warnings" dropped, all but 1Dt100h and 1Dt2d could not be given their own species

#                           #--------#
#                           |!!NOTE!!|
#                           #--------#
# I am thinking I might be able to compare msa_percent to completeness
# P.s. I had fun making that note box


combi <- gtdbtk %>% 
  rename_with(~ paste0("GTDBTK_", .x), -accession) %>% 
  left_join(
    mlst_table %>% rename_with(~ paste0("MLST2_", .x), -accession),
    by = "accession"
    ) 
# I admit I used AI help with this one, but it was too fancy and cool to not do
# It also told me that the ~ in the rename_with() statements
# is a shorthand way of doing a function, the .x is just apparently the needed formatting

rm(gtdbtk)
rm(mlst_table)
#cleaning

# 03. bring in checkm2 stuff ----------------------------------------------


checkm2 <- read_delim('../01_inputs/previous_analysis_results/checkm2/quality_report.tsv')

checkm2$accession <- gsub("flye_asm_(.*)_part2", "\\1", checkm2$Name)

combi <- combi  %>% 
  left_join(
    checkm2 %>% select(
      accession, Completeness, Contamination) %>% 
      rename_with(~ paste0("Checkm2_", .x), -accession),
    by = "accession"
  ) 
# removed columns from checkm2: completeness_model_used, translation_table, additional_notes
# model was same for all but 100h
rm(checkm2)
# 04. NCBI REST API -------------------------------------------------------


ncbi_accessions <- c(combi$GTDBTK_closest_placement_reference, combi$MLST2_Reference_assembly)
ncbi_accessions <- ncbi_accessions[ncbi_accessions != "N/A"]
ncbi_accessions <- gsub(".1$", "", ncbi_accessions)
ncbi_accessions <- unique(ncbi_accessions)
ncbi_accessions
# filtering and cleaning the columns so that I have a neat list to pass to the API

urlstring <- 'https://api.ncbi.nlm.nih.gov/datasets/v2/'
api_header <- read.delim("../98_config/privatestuff.csv", header = T, sep = ",") %>%
  filter(username == "ncbi") %>% pull(password) %>% as.character()
# API key is private, so need to bring it in from place not on public repo
#                           #--------#
#                           |!!NOTE!!|
#                           #--------#
# needs improvement, had to get AI help, not the best it could be

ncbi_accessions <- paste0(ncbi_accessions, collapse = ",")
# needed for API to work

req <- request(urlstring) %>% 
  req_auth_bearer_token(api_header) %>% 
  req_url_path_append('genome', 'accession', ncbi_accessions, 'dataset_report') %>% 
  req_headers(Accept = 'application/json') %>% 
  req_perform(path = paste0("../03_outputs/autoMLST2_subproject/metadata/reference_accessions_metadata_", "{Sys.Date()}.json"))
# API call for dataset_report to give metadata, puts it in a landing pad, {Sys.Date()} should add date, but doesnt work
# formatting got from httr2 documentation

# read file back in so req doesnt have to be run more than once
json_data <- jsonlite::read_json("../03_outputs/autoMLST2_subproject/metadata/reference_accessions_metadata_2025-05-27.json")

reports_list <- json_data$reports %>% 
  enframe() %>% 
  unnest_wider(value) %>% 
  hoist(assembly_info,
        "assembly_status",
        "refseq_category",
        "biosample") %>%  
  hoist(biosample,
        "geo_loc_name",
        "lat_lon",
        "host",
        "isolation_source") %>% 
  select(accession, assembly_status, refseq_category, geo_loc_name, host, isolation_source)
# pulls out important stuff because I get a monster amount of data back
# then filters to jsut the stuff Aaron asked for
# once again we see that several groups forget to put the info in
# p.s. one of the samples has a host of harbour seal, love it

# before I can add reports_list to combi, i need to modify it slightly
reports_list$alias <- gsub(".1$", "", reports_list$accession)
# create a column without .1 end for easier merge
reports_list <- reports_list %>%   unite("host_isolation_source", host, isolation_source, sep = " / ", na.rm = TRUE, remove = FALSE)
# this combines the host and isolation source into one column for clarity

# NEXT STEP: add this into combi, in both place 1(gtdbtk) and place2 (mlst2)

combi <- combi  %>% 
  left_join(
    reports_list %>% select(
      accession, assembly_status, refseq_category, geo_loc_name, host_isolation_source) %>% 
      rename_with(~ paste0("GTDBTK_", .x), -accession),
    by = join_by(GTDBTK_closest_placement_reference == accession)
  ) %>% 
  left_join(
    reports_list %>% select(
      alias, assembly_status, refseq_category, geo_loc_name, host_isolation_source) %>% 
      rename_with(~ paste0("MLST2_", .x), -alias),
    by = join_by(MLST2_Reference_assembly == alias)
  ) %>% 
  select(
    accession, 
    starts_with("GTDBTK_"),
    starts_with("MLST2_"),
    everything()
  )
# this adds on the data twice for the two methods, there is some overlap
# needed some AI help for that final select, i was close tho
# I thought there was a function called "begins_with", so I was close

rm(api_header, ncbi_accessions, urlstring, json_data, reports_list, req)
# cleaning

# 05. Tabulation Beautification -------------------------------------------

# before the table, I need to do some display names, because doing 19 columns
# of names manually would be painful
combi_names <- data.frame(orig_name = names(combi), new_name = 
                            c("Accession", "Closest Taxonomy Placement", 
                              "Closest Placement Reference",
                "Closest Placement ANI", "MSA%", "Assembly Status", "Reference Genome",
                "Location Name", "Host / Isolation Source", "Reference Assembly",
                "Reference Name", "MASH Distance", "Estimated ANI", "Assembly Status",
                "Reference Genome", "Location Name", "Host / Isolation Source",
                "Completeness", "Contamination"))

# needs to be a vector according to documentation
rename_vector <- setNames(combi_names$new_name, combi_names$orig_name)
# now that I have the names, I just need to update the frame before tabling

# need to replace "N/A" with " - " for clarity
combi <- combi %>% 
  mutate(across(where(is.character), ~ replace(.x, .x == "N/A", NA)))
# used AI for this, it turns "N/A" into a true NA, so it will be picked up later.

# got the table code from previous use in a prac write up

frog_flextable <- flextable(combi) %>% 
  # general
    theme_vanilla() %>%
  # header
  set_header_labels(values = rename_vector) %>%
    add_header_row(colwidths = c(1, 8, 8, 2),
                   values = c("", "GTDB-Tk", "autoMLST2", "CheckM2")) %>% 
    align(align = "center", 
          part = "header") %>%
  # body 
  compose(j = "GTDBTK_refseq_category", 
          value = mk_par(
            ifelse(frog_flextable$body$dataset$GTDBTK_refseq_category == "reference genome", "Y", "")
          )) %>% 
  colformat_char(na_str = " - ") %>% 
  # NEXT STEPS: also, colour numbers in / replace reference genome with Y/N
    #footer stuff
  vline(i = NULL, j = 1, border = fp_border_default(), part = "all") %>% 
  vline(i = NULL, j = 9, border = fp_border_default(), part = "all") %>%
  vline(i = NULL, j = 17, border = fp_border_default(), part = "all") %>%
  italic(j = c("GTDBTK_closest_tax", "MLST2_Ref_name")) %>% 
  # footer
    add_footer_lines("Data collected on 2025-03-10 at Deiniol Road, Brambell Building, 1st Floor Lab B1") %>% 
    color(part = "footer", color = "#666666") %>%
    
    footnote(i = 1, j = 2,
             part = "header",
             ref_symbols = "*",
             value = as_paragraph("P values taken between pairs of times when measurements were taken for a Wilcoxon Paired Rank Test")) %>% 
    # Apply to all numeric columns
    color(
      color = color_if_greater,
      part = "body",
      j = numeric_cols
    )

frog_flextable
combi_names "MLST2_refseq_category"
combi$GTDBTK_refseq_category
