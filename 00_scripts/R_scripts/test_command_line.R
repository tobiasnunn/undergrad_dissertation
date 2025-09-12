#simple R script to load RDS file and print first 3 list objects
# Usage: Rscript simple_rds_reader.R filename.rds

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if filename was provided
if (length(args) == 0) {
  stop("Please provide an RDS filename as an argument")
}

filename <- args[1]

# Check if file exists
if (!file.exists(filename)) {
  stop(paste("File", filename, "does not exist"))
}

# Load the RDS file
cat("Loading RDS file:", filename, "\n")
data_list <- readRDS(filename)

# Check if it's a list
if (!is.list(data_list)) {
  stop("The RDS file does not contain a list object")
}

# Print information about the list
cat("List contains", length(data_list), "objects\n")
cat("Names of objects:", paste(names(data_list), collapse = ", "), "\n\n")

# Print head of first 3 objects
n_objects <- min(3, length(data_list))

for (i in 1:n_objects) {
  cat("=== Object", i, "===\n")
  if (!is.null(names(data_list)[i])) {
    cat("Name:", names(data_list)[i], "\n")
  }
  cat("Class:", class(data_list[[i]])[1], "\n")
  print(head(data_list[[i]]))
  cat("\n")
}
