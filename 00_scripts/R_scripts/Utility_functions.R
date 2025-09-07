# these are functions for processes I perform in multiple scripts


# 01. load private information --------------------------

get_secrets <- function(secret_name){
  config_file <- "../98_config/privatestuff_v2.csv" # relative path to file
  
  if (!file.exists(config_file)) { # if file does not exist, do error handling
    stop(paste("Config file", config_file, "not found"))
  }
  
  # read in the secrets file
  secrets_file <- read.delim("../98_config/privatestuff_v2.csv", header = T, sep = ",")
  # filter the file so I either get the postgres secrets or the NCBI secrets (this is secret_name)

  # more error handling, incase an invalid secret_name is supplied, i.e. "banana"
  if (!secret_name %in% secrets_file$purpose) { 
    stop(paste("Secret name", secret_name, "not found"))
  } 
  
  # find match in file with secret_name
  secret_row <- secrets_file[which(secrets_file$purpose == secret_name), ]
  secret_list <- as.list(secret_row[, !names(secret_row) %in% "purpose"])
  
  return(secret_list)
}


# 02. check if file path exists (if not, create it)-------------------------------------------

check_and_create_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    message("Directory ", dir_path, " does not exist. Creating...")
    
    # Try to create directory
    created <- dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    
    if (!created) {
      stop("Failed to create directory: ", dir_path)
    }
    
    message("Directory created successfully: ", dir_path)
  }
}
