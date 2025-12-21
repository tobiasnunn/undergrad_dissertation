# description: Compacts DuckDB database by copying to new file and removing unused space

library(duckdb)
old_db <- "02_analysis/analysis_database/pseudomonas.db"
new_db <- "02_analysis/analysis_database/new.db"
backup_db <- "02_analysis/analysis_database/existing.db.backup"

# Remove new.db if it exists
if (file.exists(here::here(new_db))) file.remove(here::here(new_db))

# Connect to new database
con <- dbConnect(duckdb(), here::here(new_db))

# Attach old database and copy everything
dbExecute(con, sprintf("ATTACH '%s' AS existing", old_db))
dbExecute(con, "COPY FROM DATABASE existing TO new")
dbExecute(con, "DETACH existing")

dbDisconnect(con, shutdown = TRUE)

# Replace old with new
file.rename(here::here(old_db), here::here(backup_db))
file.rename(here::here(new_db), here::here(old_db))

# Clean up WAL files
wal_files <- c(paste0(here::here(old_db), ".wal"), 
               paste0(here::here(old_db), ".wal.tmp"),
               paste0(here::here(backup_db), ".wal"),
               paste0(here::here(backup_db), ".wal.tmp"))

for (wal in wal_files) {
  if (file.exists(wal)) file.remove(wal)
}

cat("Database compacted!\n")
file.remove(here::here(backup_db))
