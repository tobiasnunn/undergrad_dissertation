
# library -----------------------------------------------------------------

library(flextable)
library(duckdb)
library(tidyverse)


# bring in tables from duckdb ---------------------------------------------

# load stuff needed for enrichKO
con <- dbConnect(duckdb(), here::here("02_analysis/analysis_database/pseudomonas.db"), read_only=TRUE)

# get summary table for filtering steps
filter_table <- dbGetQuery(con, 
                                  "SELECT * FROM animal_hosts;")


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



# table of host counts
host_group_counts <- dbGetQuery(con, 
"SELECT animal_group as 'Animal group', COUNT(*) as Count
FROM (SELECT DISTINCT Accession, animal_group
FROM pseudomonas_annotations pa) as Q1
GROUP BY animal_group;")

#table of taxonomy HERE, maybe best to do both, a host group broken down by taxonomy


# close connection to duckdb
dbDisconnect(con)


# Supplementary Table 1 - filter ------------------------------------------

# keeping the non-used groups bloats the table, maybe useful, but doesnt do anything for the analysis
filter_table_filtered <- filter(filter_table, filter_table$animal_group %in% 
                                  c("Amphibian", "Bird", "Fish", "Mammal")) 

S1 <- flextable(filter_table_filtered) %>% 
  set_caption("Supplementary Table 1: Step to filter out Domesticated hosts") %>% 
  add_footer_lines("Step taken to remove domesticated hosts, as the microbiomes they harbor may not be representative of wild hosts.")

S1

# will need to do more tables, 1 for the Country filter I had Claude do and 1 for the iso source filter I am about to have him do.



# Supplementary Table 2 - Precise host distribution -----------------------

s2 <- flextable(host_group_counts) %>% 
  set_caption("Supplementary Table 2: Precise distribution of host grouping amongst study samples")

s2

# Supplementary Table 3 - Precise taxonomic info --------------------------


# S table 4 - study pseudomonas -------------------------------------------

s4 <- flextable(study_pseudomonas) %>% 
  set_caption("Supplementary Table 4: Full key information relating to study samples")

s4
# save_as_image(s4, "03_outputs/03_descriptive_tables/study_pseudomonas.png")
# maybe better to save as tsv

# make the isolation source information -----------------------------------


iso_source <- select(study_pseudomonas, animal_group, isolationSource)
iso_source
# send to Claude

#code also from Claude
iso_source <- iso_source %>%
  mutate(
    sample_location = case_when(
      # INTERNAL - Respiratory
      grepl("lung|trachea|respiratory|bronchial|BAL", isolationSource, ignore.case = TRUE) ~ "internal_respiratory",
      
      # INTERNAL - Urogenital
      grepl("uterus|urinary|urine|genital|vagina|clitoris|Blowhole", isolationSource, ignore.case = TRUE) ~ "internal_urogenital",
      
      # INTERNAL - Gastrointestinal
      grepl("intestine|gut|rectum|jejunum|feces|fecal|gastrointestinal", isolationSource, ignore.case = TRUE) ~ "internal_gi",
      
      # INTERNAL - Other organs/fluids
      grepl("spleen|organ|blood|milk|semen|saliva", isolationSource, ignore.case = TRUE) ~ "internal_other",
      
      # INTERNAL - Gills (fish-specific respiratory)
      grepl("gill", isolationSource, ignore.case = TRUE) ~ "internal_gills",
      
      # EXTERNAL - Ear
      grepl("ear|oto", isolationSource, ignore.case = TRUE) ~ "external_ear",
      
      # EXTERNAL - Skin/wounds
      grepl("skin|wound|lesion|ulcer|nail pus|elbow", isolationSource, ignore.case = TRUE) ~ "external_skin",
      
      # EXTERNAL - Cloaca
      grepl("cloaca", isolationSource, ignore.case = TRUE) ~ "external_cloaca",
      
      # EXTERNAL - Eye
      grepl("eye|ocular", isolationSource, ignore.case = TRUE) ~ "external_eye",
      
      # EXTERNAL - Other appendages/surfaces
      grepl("wing|foot", isolationSource, ignore.case = TRUE) ~ "external_appendage",
      
      # EXTERNAL - Swabs (rectal, skin swabs not caught above)
      grepl("swab", isolationSource, ignore.case = TRUE) ~ "external_swab",
      
      # UNKNOWN - Environmental/context
      grepl("farm|facility|lake|necropsy|collected|applicable", isolationSource, ignore.case = TRUE) ~ "unknown_environmental",
      
      # UNKNOWN - Missing
      is.na(isolationSource) | isolationSource == "missing" ~ "unknown_missing",
      
      # UNKNOWN - Ambiguous terms
      isolationSource %in% c("fish", "alevin") ~ "unknown_ambiguous",
      
      # Catch-all
      TRUE ~ "unknown_other"
    )
  )

# Now you can see both columns side by side
iso_source %>%
  select(animal_group, isolationSource, sample_location) %>%
  head(20)
unique(iso_source$sample_location)
