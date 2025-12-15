
# 01_packages -------------------------------------------------------------


# packages for doing database connections
library(RPostgres)
library(DBI)


# 02_establish_connection -------------------------------------------------


# reads in my private credentials for my local Postgres server
# this file is not uploaded to the repo for obvious reasons
private <- read.csv("../98_config/privatestuff.csv", sep=",", header=T)

# Connect to a specific postgres database
conn <- dbConnect(RPostgres::Postgres(), 
                  dbname = 'postgres',
                  host = 'localhost', 
                  port = 5432, 
                  user = private$username,
                  password = private$password)

# 03_create_or_connect_to_database ----------------------------------------


# Check if database exists
db_obj <- dbGetQuery(conn, "SELECT 1 FROM pg_database WHERE datname = 'genomic_analysis'")
db_exists <- nrow(db_obj) > 0

# after checking database does not exist create database named "genomic_analysis"
if (!db_exists) {
  # Database doesn't exist, create it
  # We'll use dbExecute with immediate=TRUE to execute outside a transaction
  dbExecute(conn, "CREATE DATABASE genomic_analysis", immediate = TRUE)
} else {
  cat("Database 'genomic_analysis' already exists.\n")
}

