# file to automate notebook creation
# to run, in the "terminal" section, type: "Rscript <filepath>"
# which for this file is: "../00_scripts/R_scripts/Create_notebooks.R"
# a week number must be supplied as an argument

# get the argument passed in
# trailingonly = T avoids system arguments (that I dont want) being read in.
args <- commandArgs(trailingOnly = T)

# check if arg was supplied, crash if not
if (length(args) == 0) {
  cat("no week number provided")
  quit(status = 1)
}

# retreave week number from arguments
week_number <- args[1]

# check to see if file already exists
filename <- paste0("posts/week_", week_number, ".qmd")

if (file.exists(filename)) {
  cat(filename, " already exists")
  quit(status = 1)
}

#check the current date so that the publish date is correct
# (I always set it to the Sunday of the week)
next_sunday <- as.character(Sys.Date() + (7 - as.numeric(format(Sys.Date(), "%u"))))
last_monday <- as.character(Sys.Date() - (as.numeric(format(Sys.Date(), "%u")) - 1))
# create content of new notebook
content <- paste0(
"---
title: \"Week ", week_number, " ", last_monday, " to ", next_sunday, "\"
author: \"", "Tobias Nunn", "\"
date: \"", next_sunday, "\"
categories:
---

# Introduction

*Brief overview of the week's objectives and context.*


# Methods

*Description of stuff I did.*


# Results

*Key findings, observations, and data from this week's work.*


# Conclusion

*Summary of outcomes, implications, and next steps.*


# Scientific Papers

*Relevant literature read this week.*

-
-
-

")

#write the content to a new file

writeLines(content, filename)
cat("File created")