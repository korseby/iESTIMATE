


# ############################## iESTIMATE ##############################



# ---------- Preparations ----------
library(devtools)
library(roxygen2)
library(usethis)
library(testthat)



# ---------- Create new R package ----------
setwd("./")

# Create new package
#usethis::create_package("iESTIMATE")
#devtools::create(path="./")

# Add function from R files
usethis::use_r("functions.R")

# Load functions
load_all()

# Check requirements
devtools::check()



