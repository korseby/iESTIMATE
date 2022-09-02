#!/usr/bin/env Rscript



# ---------- Preparations ----------
# Load libraries
library(parallel)               # Detect number of cpu cores
library(foreach)                # For multicore parallel
library(doMC)                   # For multicore parallel
library(RColorBrewer)           # For colors
library(multtest)               # For diffreport
library(MSnbase)                # MS features
library(xcms)                   # Swiss army knife for metabolomics
library(CAMERA)                 # Metabolite Profile Annotation
library(Spectra)                # Spectra package needed for XCMS3
library(mixOmics)               # Statistics for various Omics
# ERROR: vegan conflicts with mda and klaR, unload packages before using any of the analyses !!!
if ("package:mda" %in% search()) detach(package:mda, unload=TRUE)
if ("package:klaR" %in% search()) detach(package:klaR, unload=TRUE)
library(vegan)
library(multcomp)               # For Tukey test
library(Hmisc)                  # For correlation test
library(gplots)                 # For fancy heatmaps
library(lme4)                   # Linear Models with Random Effects
library(lmerTest)               # Create p-values for LMER
library(ape)                    # Phylogeny
library(pvclust)                # Phylogeny
library(dendextend)             # Phylogeny
library(cba)                    # Phylogeny
library(phangorn)               # Phylogeny
library(ontologyIndex)          # Reading obo ontology files
library(webchem)                # Converting InChI to InChIKey
library(circlize)               # For sunburst plot
library(plotrix)                # For sunburst plot
library(Boruta)                 # Random-Forest BORUTA
library(rpart)                  # Regression trees
library(caret)                  # Swiss-army knife for statistics
library(pROC)                   # Evaluation metrics
library(PRROC)                  # Evaluation metrics
library(multiROC)               # Evaluation metrics

# Setup R error handling to go to stderr
#options(show.error.messages=F, error=function() { cat(geterrmessage(), file=stderr()); q("no",1,F) } )

# Set locales and encoding
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
loc <- Sys.setlocale(category="LC_ALL", locale="C")
options(encoding="UTF-8")

# Set options
options(stringAsfactors=FALSE, useFancyQuotes=FALSE)

# Multicore parallel
#nSlaves <- detectCores(all.tests=FALSE, logical=FALSE)
#nSlaves <- 16
#registerDoMC(nSlaves)



# ---------- User variables ----------
# Take in trailing command line arguments
#args <- commandArgs(trailingOnly=TRUE)
#if (length(args) < 4) {
#    print("Error! No or not enough arguments given.")
#    print("Usage: $0 metfrag_results.csv metfrag_hits_limit metfrag_structure_method output.pdf")
#    quit(save="no", status=1, runLast=FALSE)
#}



# Working directory
setwd("~/Desktop/Projekte/Habilitation/Mosses/iestimate_mtbls709/")

# Data directory
mzml_dir <- "~/Desktop/Projekte/Habilitation/Mosses/ms2-dda/data/raw/"

# MS1 variables
polarity <- "positive"
pol <- substr(x=polarity, start=1, stop=3)
ppm <- 35
peakwidth <- c(4,21)
prefilter <- c(5,50)
#fitgauss <- FALSE
#verbose.columns <- FALSE
snthresh <- 11
mzdiff <- 0.01
noise <- 0
integrate <- 1
minfrac <- 0.7
bwindow <- 0.25
minfracspan <- 0.2
ms1_intensity_cutoff <- 100	          # approx. 0.01% (= log 100)

# Peak detection variables

# MS2 variables
mzabs <- 0.01                         # Absolute mass error (in seconds) used for merging MS/MS spectra
mzppm <- 5 #5                         # ppm error used for merging MS/MS spectra
rtabs <- 5 #20                        # Retention time error (in seconds) used for merging MS/MS spectra
max.rt.range <- 20                    # Permitted retention time window (in seconds) of grouped MS1 precursors
max.mz.range <- 0.01                  # Permitted m/z window of grouped MS1 precursors
min.rt <- 10                          # Minimum retention time for selected precursors
max.rt <- 1020                        # Maximum retention time for selected precursors
min.mz <- 50                          # Minimum m/z value for selected precursors
max.mz <- 1500                        # Maximum m/z value for selected precursors
msms.intensity.threshold <- 100       # Minimum intensity value for MS/MS peaks

# Classifier options
mzabs <- 0.01                        # Absolute mass error (in seconds) used for merging MS/MS spectra
mzppm <- 5 #35                       # ppm error used for merging MS/MS spectra
msms.intensity.threshold <- 100 #10  # Minimum intensity value for MS/MS peaks

CANOPUS <- TRUE                      # Use CANOPUS classification instead of MetFamily
SIRIUS_VERSION <- 5

# Preparations for plotting
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)

# Save and load RData
#save.image()
load(".RData")



# ---------- Peak detection ----------
source("_functions.r")

setMSnbaseFastLoad(TRUE)

source("peak_detection_pos.r")



