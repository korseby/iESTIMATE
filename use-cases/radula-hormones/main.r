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
setwd("C:/Users/kaitl/Desktop/Research/UNB/Research Project/Germany/Bioinformatics/Radula_hormones_CANOPUS")

# Data directory
mzml_dir <- "C:/Users/kaitl/Desktop/Research/UNB/Research Project/Germany/Bioinformatics/Radula_hormones_CANOPUS/raw/"

# MS1 variables
polarity <- "positive"
pol <- substr(x=polarity, start=1, stop=3)
ppm <- 35
ms1_intensity_cutoff <- 100	          # approx. 0.01% (= log 100)

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
mzppm <- 35                          # ppm error used for merging MS/MS spectra
msms.intensity.threshold <- 10       # Minimum intensity value for MS/MS peaks

CANOPUS <- TRUE                      # Use CANOPUS classification instead of MetFamily

# Preparations for plotting
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)

# Save and load RData
#save.image()
#load(".RData")



# ---------- Peak detection ----------
source("_functions-1.r")

setMSnbaseFastLoad(TRUE)

source("peak_detection_pos.r")
source("peak_detection_neg.r")



# ---------- Create merged pos+neg objects ----------
# factors
mzml_pheno_samples <- mzml_pheno_samples_pos
mzml_pheno_colors <- mzml_pheno_colors_pos

# feat_list
feat_list <- cbind(feat_list_pos, feat_list_neg)
rownames(feat_list) <- gsub(x=rownames(feat_list_pos), pattern="\\.auto.*", replacement="")
colnames(feat_list) <- c(paste0(colnames(feat_list_pos),"_pos"), paste0(colnames(feat_list_neg),"_neg"))

# bina_list
bina_list <- cbind(bina_list_pos, bina_list_neg)
rownames(bina_list) <- gsub(x=rownames(bina_list_pos), pattern="\\.auto.*", replacement="")
colnames(bina_list) <- c(paste0(colnames(bina_list_pos),"_pos"), paste0(colnames(bina_list_neg),"_neg"))

# uniq_list
uniq_list <- cbind(uniq_list_pos, uniq_list_neg)
rownames(uniq_list) <- gsub(x=rownames(uniq_list_pos), pattern="\\.auto.*", replacement="")
colnames(uniq_list) <- c(paste0(colnames(uniq_list_pos),"_pos"), paste0(colnames(uniq_list_neg),"_neg"))

# Create data frame
model_div             <- data.frame(features=apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div$richness    <- apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } )
#model_div$menhinick   <- apply(X=bina_list, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div$shannon     <- apply(X=feat_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div$pielou      <- apply(X=scale(feat_list, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list, species), index="chao")
model_div$simpson     <- apply(X=feat_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div$inverse     <- apply(X=feat_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div$fisher      <- apply(X=feat_list, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div$unique      <- apply(X=uniq_list, MARGIN=1, FUN=function(x) { sum(x) })

# Remove NAs if present
model_div[is.na(model_div)] <- 0

# class_list
class_list <- cbind(class_list_pos, class_list_neg)
rownames(class_list) <- gsub(x=rownames(class_list_pos), pattern="\\.auto.*", replacement="")
class_list <- cbind(sapply(unique(colnames(class_list)[duplicated(colnames(class_list))]), function(x) rowSums(class_list[,grepl(paste(x, "$", sep=""), colnames(class_list))])), class_list[,!duplicated(colnames(class_list)) & ! duplicated(colnames(class_list), fromLast=TRUE)])

# superclass_list
superclass_list <- cbind(superclass_list_pos, superclass_list_neg)
rownames(superclass_list) <- gsub(x=rownames(superclass_list_pos), pattern="\\.auto.*", replacement="")
superclass_list <- cbind(sapply(unique(colnames(superclass_list)[duplicated(colnames(superclass_list))]), function(x) rowSums(superclass_list[,grepl(paste(x, "$", sep=""), colnames(superclass_list))])), superclass_list[,!duplicated(colnames(superclass_list)) & ! duplicated(colnames(superclass_list), fromLast=TRUE)])

# Create a div_classes and div_superclasses for sunburst
onlyPos <- row.names(div_classes_pos)[!row.names(div_classes_pos)%in% row.names(div_classes_neg)] # only in 'first', same as: setdiff(first, second)
onlyNeg <- row.names(div_classes_neg)[!row.names(div_classes_neg) %in% row.names(div_classes_pos)] # only in 'second', same as: setdiff(second, first)


onlyPos <- data.frame(
	row.names = onlyPos,
	frequency = rep(0, length(onlyPos))
)
onlyNeg <- data.frame(
	row.names = onlyNeg,
	frequency = rep(0, length(onlyNeg))
)

div_classes_pos_prep <- rbind(div_classes_pos, onlyNeg)
div_classes_neg_prep <- rbind(div_classes_neg, onlyPos)
div_classes <- merge(div_classes_pos_prep, div_classes_neg_prep, by='row.names', all=TRUE)
div_classes <- data.frame(
	row.names = div_classes$Row.names,
	frequency.x = div_classes$frequency.x, 
	frequency.y = div_classes$frequency.y
)
div_classes <- rowSums(div_classes)
div_classes <- as.data.frame(div_classes)
colnames(div_classes) <- "frequency"

# superclass level
library(tidyr)
div_superclasses <- data.frame(class = row.names(div_classes), frequency = div_classes$frequency)
div_superclasses <-	separate(div_superclasses, col = class, into = c("Organic", "Superclass", "Class"), sep = ";")
library(dplyr)
div_superclasses <- div_superclasses %>% 
	group_by(Superclass) %>% 
	summarise(frequency = sum(frequency))
div_superclasses <- as.data.frame(div_superclasses)
row.names(div_superclasses) <- div_superclasses$Superclass

# Write a big classifier dataframe (can be used to extract boruta results)
classifier_canopus_pos$Metabolite.name <- paste0(classifier_canopus_pos$Metabolite.name,"_pos")
classifier_canopus_neg$Metabolite.name <- paste0(classifier_canopus_neg$Metabolite.name,"_neg")
classifier_canopus <- rbind(classifier_canopus_neg, classifier_canopus_pos)


# ############################## MS1 statistics ##############################



# ---------- Histogram ----------
pdf(file="plots/ms1_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(feat_list))
dev.off()



# ---------- Variation partitioning ----------
mzml_pheno_growth_stress_samples <- as.factor(c("Stress", "Stress", "Stress", "Stress", "Stress", "Stress", 
												"Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth",
												"Stress", "Stress", "Stress", "Stress", "Stress", "Stress",
												"Growth", "Growth", "Growth", "Growth", "Growth", "Growth",
												"Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress"))
model_varpart <- varpart(scale(feat_list), ~ mzml_pheno_samples, ~ mzml_pheno_growth_stress_samples)

mzml_pheno_methanol_samples <- as.factor(c("1", "1", "1", "10", "10", "10", 
																					 "1", "1", "1", "10", "10", "10",
																					 "0", "0", "0",
																					 "1", "1", "1", "10", "10", "10", 
																					 "1", "1", "1", "10", "10", "10",
																					 "0", "0", "0",
																					 "1", "1", "1", "10", "10", "10"))



model_varpart <- varpart(scale(feat_list), ~ mzml_pheno_samples, ~ mzml_pheno_methanol_samples, ~mzml_pheno_growth_stress_samples)

# Plot results
pdf(file="plots/ms1_varpart_3variables.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("Samples","MeOH", "Growth Time"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","red", "purple"))
dev.off()

# PCA of combined feature table
pdf(file="plots/ms1_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca <- prcomp(feat_list, center=TRUE)
plot(ms1_pca$x[, 1], ms1_pca$x[,2], pch=19, main="PCA: Grouping of samples",
		 xlab=paste0("PC1: ", format(summary(ms1_pca)$importance[2, 1] * 100, digits=3), " % variance"),
		 ylab=paste0("PC2: ", format(summary(ms1_pca)$importance[2, 2] * 100, digits=3), " % variance"),
		 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca$x[, 1], ms1_pca$x[,2], labels=names(feat_list[, 5]), col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()


# ---------- Diversity ----------
# Plot unique features
pdf(paste("plots/ms1_richness_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$unique ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Number of unique features", xlab="treatment", ylab="number of unique features")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$unique, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot total features
pdf(paste("plots/ms1_richness_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$features ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Number of features", xlab="treatment", ylab="number of features")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$features, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Shannon index
pdf(paste("plots/ms1_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$shannon ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Pielou evenness for species
pdf(paste("plots/ms1_diversity_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$pielou ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Pielou\'s evenness", xlab="treatment", ylab="Pielou diversity index (J)")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$pielou, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()


# ---------- Variable Selection ----------
# BORUTA
#library(doParallel)
#cl <- makeCluster(16)
#registerDoParallel(cl)
sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=feat_list, sel_factor=mzml_pheno_growth_stress_samples, sel_colors=mzml_pheno_colors$cols, sel_method="best", tune_length=10, pvalue=0.01, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms1_select_boruta_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
f.heatmap.selected_features(feat_list=feat_list, sel_feat=sel_boruta_rf, filename="plots/ms1_select_boruta.pdf", main="BORUTA")
#stopCluster(cl)

#Search for MS1 spectra that have useable MS2 spectra
ms2_names_pos <- paste(names(ms2_spectra_pos), "_pos", sep = "")
ms2_names_neg <- paste(names(ms2_spectra_neg), "_neg", sep = "")
ms2_combined_names <- c(ms2_names_neg, ms2_names_pos)
mod_feat_list <- feat_list[, ms2_combined_names]

mzml_pheno_treatments_concentrations_samples <- as.factor(c("AA1", "AA1", "AA1", "AA10", "AA10", "AA10", 
																														"BAP1", "BAP1", "BAP1", "BAP10", "BAP10", "BAP10",
																													  "G", "G", "G",
																												 	  "JA1", "JA1", "JA1", "JA10", "JA10", "JA10",
																														"NAA1", "NAA1", "NAA1", "NAA10", "NAA10", "NAA10",
																														"SA1", "SA1", "SA1", "SA10", "SA10", "SA10",
																														"S", "S", "S"))

sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=mod_feat_list, sel_factor=mzml_pheno_treatments_concentrations_samples, sel_colors=mzml_pheno_colors$cols, sel_method="best", tune_length=10, pvalue=0.01, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms1_select_boruta_roc_modMS1_test1.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
f.heatmap.selected_features(feat_list=mod_feat_list, sel_feat=sel_boruta_rf, filename="plots/ms1_mod_select_boruta_test1.pdf", main="BORUTA")

# Get BORUTA class information (copy important features from PDF)
boruta_classes <- classifier_canopus[classifier_canopus$Metabolite.name %in% MS1.MOD.BORUTA.Features$V1,]
write.csv(boruta_classes, file = "data/ms1_mod_boruta_features.csv")

# Break feature table up by hormone type
stress_names <- c("AA1.1", "AA1.2", "AA1.3", "AA10.1", "AA10.2", "AA10.3",
									"JA1.1", "JA1.2", "JA1.3", "JA10.1", "JA10.2", "JA10.3",
									"SA1.1", "SA1.2", "SA1.3", "SA10.1", "SA10.2", "SA10.3",
									"S.Control.1", "S.Control.2", "S.Control.3")
mod_feat_list_stress <- mod_feat_list[row.names(mod_feat_list) %in% stress_names,]
growth_names <- c("BAP1.1", "BAP1.2", "BAP1.3", "BAP10.1", "BAP10.2", "BAP10.3",
									"NAA1.1", "NAA1.2", "NAA1.3", "NAA10.1", "NAA10.2", "NAA10.3",
									"G.Control.1", "G.Control.2", "G.Control.3")
mod_feat_list_growth <- mod_feat_list[row.names(mod_feat_list) %in% growth_names,]

stress_factor <- as.factor(gsub(stress_names, pattern = "\\..", replacement = ""))
growth_factor <- as.factor(gsub(growth_names, pattern = "\\..", replacement = ""))

feat_list_stress <- feat_list[row.names(feat_list) %in% stress_names,]
feat_list_growth <- feat_list[row.names(feat_list) %in% growth_names,]

hormone_names <- c("AA1.1", "AA1.2", "AA1.3", "AA10.1", "AA10.2", "AA10.3",
									 "BAP1.1", "BAP1.2", "BAP1.3", "BAP10.1", "BAP10.2", "BAP10.3",
									 "G.Control.1", "G.Control.2", "G.Control.3",
									"JA1.1", "JA1.2", "JA1.3", "JA10.1", "JA10.2", "JA10.3",
									"NAA1.1", "NAA1.2", "NAA1.3", "NAA10.1", "NAA10.2", "NAA10.3",
									"SA1.1", "SA1.2", "SA1.3", "SA10.1", "SA10.2", "SA10.3",
									"S.Control.1", "S.Control.2", "S.Control.3")

hormone_factors <- as.factor(gsub(hormone_names, pattern = "\\..", replacement = ""))


# Boruta for stress
sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=feat_list_stress, sel_factor = stress_factor, sel_colors=mzml_pheno_colors$cols, sel_method="best", tune_length=10, pvalue=0.05, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms1_select_boruta_roc_MS1_stress.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
print(sel_boruta_rf[['_RMSE_']])
print(sel_boruta_rf[['_R2_']])
print(sel_boruta_rf[['_MAE_']])
f.heatmap.selected_features(feat_list=mod_feat_list_stress, sel_feat=sel_boruta_rf, filename="plots/ms1_select_boruta_stress.pdf", main="BORUTA")

# Boruta for growth
sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=feat_list_growth, sel_factor=growth_factor, sel_colors=mzml_pheno_colors$cols, sel_method="best", tune_length=2, pvalue=0.05, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms1_select_boruta_roc_MS1_growth.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
print(sel_boruta_rf[['_RMSE_']])
print(sel_boruta_rf[['_R2_']])
print(sel_boruta_rf[['_MAE_']])
f.heatmap.selected_features(feat_list=mod_feat_list_growth, sel_feat=sel_boruta_rf, filename="plots/ms1_select_boruta_growth.pdf", main="BORUTA")


# Normalize Feature Table
stress_treatments <- c("AA1.1", "AA1.2", "AA1.3", "AA10.1", "AA10.2", "AA10.3",
											 "JA1.1", "JA1.2", "JA1.3", "JA10.1", "JA10.2", "JA10.3",
										   "SA1.1", "SA1.2", "SA1.3", "SA10.1", "SA10.2", "SA10.3")
growth_treatments <- c("BAP1.1", "BAP1.2", "BAP1.3", "BAP10.1", "BAP10.2", "BAP10.3",
											 "NAA1.1", "NAA1.2", "NAA1.3", "NAA10.1", "NAA10.2", "NAA10.3")

feat_list_stress_treatments <- feat_list[row.names(feat_list) %in% stress_treatments,]
feat_list_growth_treatments <- feat_list[row.names(feat_list) %in% growth_treatments,]

stress_controls <- c("S.Control.1", "S.Control.2", "S.Control.3")
growth_controls <- c("G.Control.1", "G.Control.2", "G.Control.3")

feat_list_stress_controls <- feat_list[row.names(feat_list) %in% stress_controls,]
feat_list_growth_controls <- feat_list[row.names(feat_list) %in% growth_controls,]

feat_list_stress_controls_mean <- colMeans(as.data.frame(feat_list_stress_controls))
feat_list_growth_controls_mean <- colMeans(as.data.frame(feat_list_growth_controls))

feat_list_stress_treatments_normalized <- feat_list_stress_treatments
i <- 1
for(val in stress_treatments){
	feat_list_stress_treatments_normalized[i, ] = (feat_list_stress_treatments[i, ] - feat_list_stress_controls_mean)/feat_list_stress_controls_mean
	i = i+1
}

feat_list_growth_treatments_normalized <- feat_list_growth_treatments
i <- 1
for(val in growth_treatments){
	feat_list_growth_treatments_normalized[i, ] = (feat_list_growth_treatments[i, ] - feat_list_growth_controls_mean)/feat_list_growth_controls_mean
	i = i+1
}

feat_list_normalized <- rbind(feat_list_growth_treatments_normalized, feat_list_stress_treatments_normalized)

hormone_names_normalized <- c("BAP1.1", "BAP1.2", "BAP1.3", "BAP10.1", "BAP10.2", "BAP10.3",
									 "NAA1.1", "NAA1.2", "NAA1.3", "NAA10.1", "NAA10.2", "NAA10.3",
									 "AA1.1", "AA1.2", "AA1.3", "AA10.1", "AA10.2", "AA10.3",
									 "JA1.1", "JA1.2", "JA1.3", "JA10.1", "JA10.2", "JA10.3",
									 "SA1.1", "SA1.2", "SA1.3", "SA10.1", "SA10.2", "SA10.3")

hormone_factors_normalized <- as.factor(gsub(hormone_names_normalized, pattern = "\\..", replacement = ""))
mzml_pheno_colors_normalized <- subset(mzml_pheno_colors, cols != "red" & cols != "blue")

# normalized full BORUTA
sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=feat_list_normalized, sel_factor=hormone_factors_normalized, sel_colors=mzml_pheno_colors_normalized$cols, sel_method="best", tune_length=10, pvalue=0.01, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms1_select_boruta_roc_normalized.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
print(sel_boruta_rf[['_RMSE_']])
print(sel_boruta_rf[['_R2_']])
print(sel_boruta_rf[['_MAE_']])
f.heatmap.selected_features(feat_list=feat_list_normalized, sel_feat=sel_boruta_rf, filename="plots/ms1_select_boruta_normalized.pdf", main="BORUTA")

#Search for MS1 spectra that have useable MS2 spectra
ms2_names_pos <- paste(names(ms2_spectra_pos), "_pos", sep = "")
ms2_names_neg <- paste(names(ms2_spectra_neg), "_neg", sep = "")
ms2_combined_names <- c(ms2_names_neg, ms2_names_pos)
mod_feat_list_normalized <- feat_list_normalized[, ms2_combined_names]
mod_feat_list_stress_treatments_normalized <- feat_list_stress_treatments_normalized[, ms2_combined_names]
mod_feat_list_growth_treatments_normalized <- feat_list_growth_treatments_normalized[, ms2_combined_names]

mzml_pheno_normalized_growth_stress_samples <- as.factor(c("Stress", "Stress", "Stress", "Stress", "Stress", "Stress", 
																								"Growth", "Growth", "Growth", "Growth", "Growth", "Growth",
																								"Stress", "Stress", "Stress", "Stress", "Stress", "Stress",
																								"Growth", "Growth", "Growth", "Growth", "Growth", "Growth",
																								"Stress", "Stress", "Stress", "Stress", "Stress", "Stress"))

mzml_pheno_normalized_treatments_samples <- as.factor(c("AA", "AA", "AA", "AA", "AA", "AA", 
																												"BAP", "BAP", "BAP", "BAP", "BAP", "BAP",
																												"JA", "JA", "JA", "JA", "JA", "JA",
																												"NAA", "NAA", "NAA", "NAA", "NAA", "NAA",
																												"SA", "SA", "SA", "SA", "SA", "SA"))

mzml_pheno_normalized_treatments_concentrations_samples <- as.factor(c("AA1", "AA1", "AA1", "AA10", "AA10", "AA10", 
																												"BAP1", "BAP1", "BAP1", "BAP10", "BAP10", "BAP10",
																												"JA1", "JA1", "JA1", "JA10", "JA10", "JA10",
																												"NAA1", "NAA1", "NAA1", "NAA10", "NAA10", "NAA10",
																												"SA1", "SA1", "SA1", "SA10", "SA10", "SA10"))

source("_functions-1.r")
sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=mod_feat_list_normalized, sel_factor=mzml_pheno_normalized_treatments_samples, sel_colors=mzml_pheno_colors_normalized$cols, sel_method="best", tune_length=10, pvalue=0.001, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms1_select_boruta_roc_modMS1_normalized_test23.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
print(sel_boruta_rf[['_RMSE_']])
print(sel_boruta_rf[['_R2_']])
print(sel_boruta_rf[['_MAE_']])
f.heatmap.selected_features(feat_list=mod_feat_list_normalized, sel_feat=sel_boruta_rf, filename="plots/ms1_mod_select_boruta_normalized_test23.pdf", main="BORUTA")

sel_boruta_rf$`_model_r2_`
sel_boruta_rf$`_R2_`
sel_boruta_rf$`_RMSE_`
mean(sel_boruta_rf$`_multiclass_metrics_`$MCC)

growth_norm_factor <- as.factor(gsub(growth_treatments, pattern = "\\..", replacement = ""))
stress_norm_factor <- as.factor(gsub(stress_treatments, pattern = "\\..", replacement =  ""))

#growth
sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=mod_feat_list_growth_treatments_normalized, sel_factor=growth_norm_factor, sel_colors=mzml_pheno_colors_normalized$cols, sel_method="best", tune_length=2, pvalue=0.05, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms1_select_boruta_roc_modMS1_normalized_growth.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
print(sel_boruta_rf[['_RMSE_']])
print(sel_boruta_rf[['_R2_']])
print(sel_boruta_rf[['_MAE_']])
f.heatmap.selected_features(feat_list=mod_feat_list_growth_treatments_normalized, sel_feat=sel_boruta_rf, filename="plots/ms1_mod_select_boruta_normalized_growth_test.pdf", main="BORUTA")

#stress
sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=mod_feat_list_stress_treatments_normalized, sel_factor=stress_norm_factor, sel_colors=mzml_pheno_colors_normalized$cols, sel_method="best", tune_length=2, pvalue=0.05, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms1_select_boruta_roc_modMS1_normalized_stress.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
print(sel_boruta_rf[['_RMSE_']])
print(sel_boruta_rf[['_R2_']])
print(sel_boruta_rf[['_MAE_']])
f.heatmap.selected_features(feat_list=mod_feat_list_stress_treatments_normalized, sel_feat=sel_boruta_rf, filename="plots/ms1_mod_select_boruta_normalized_stress_test.pdf", main="BORUTA")


# ############################## MS2 Classes statistics ##############################


# ---------- Histogram ----------
pdf(file="plots/ms2_classes_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(unlist(class_list)))
dev.off()



# ---------- Variation partitioning ----------
mzml_pheno_growth_stress_samples <- as.factor(c("Stress", "Stress", "Stress", "Stress", "Stress", "Stress", 
												"Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth",
												"Stress", "Stress", "Stress", "Stress", "Stress", "Stress",
												"Growth", "Growth", "Growth", "Growth", "Growth", "Growth",
												"Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress"))
model_varpart <- varpart(class_list, ~ mzml_pheno_samples, ~ mzml_pheno_growth_stress_samples)

mzml_pheno_methanol_samples <- as.factor(c("1", "1", "1", "10", "10", "10", 
																					 "1", "1", "1", "10", "10", "10",
																					 "0", "0", "0",
																					 "1", "1", "1", "10", "10", "10", 
																					 "1", "1", "1", "10", "10", "10",
																					 "0", "0", "0",
																					 "1", "1", "1", "10", "10", "10"))
model_varpart <- varpart(scale(feat_list), ~ mzml_pheno_samples, ~ mzml_pheno_methanol_samples)

# Plot results
pdf(file="plots/ms2_classes_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","Hormone"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- Diversity ----------
# Plot unique features
pdf(paste("plots/ms2_classes_richness_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$unique ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Number of unique features", xlab="treatment", ylab="number of unique features")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$unique, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot total features
pdf(paste("plots/ms2_classes_richness_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$features ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Number of features", xlab="treatment", ylab="number of features")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$features, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Shannon index
pdf(paste("plots/ms2_classes_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$shannon ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Pielou evenness for species
pdf(paste("plots/ms2_classes_diversity_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$pielou ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Pielou\'s evenness", xlab="treatment", ylab="Pielou diversity index (J)")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$pielou, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Sunburst plot
pdf(file="plots/ms2_classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_classes), div_classes$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_classes), "Number of spectra"=div_classes$frequency), file="plots/ms2_classes_sunburst.csv", row.names=FALSE)

# Make sunburst plots of classes detected in each polarity
onlyPos_classes <- data.frame(
	row.names = row.names(div_classes_pos)[!row.names(div_classes_pos)%in% row.names(div_classes_neg)],
	frequency = div_classes_pos[!row.names(div_classes_pos)%in% row.names(div_classes_neg),]
) 
onlyNeg_classes <- data.frame(
	row.names = row.names(div_classes_neg)[!row.names(div_classes_neg)%in% row.names(div_classes_pos)],
	frequency = div_classes_neg[!row.names(div_classes_neg)%in% row.names(div_classes_pos),]
) 
bothPolarities_classes <- data.frame(
	row.names = row.names(div_classes_neg)[row.names(div_classes_neg)%in% row.names(div_classes_pos)],
	frequency = div_classes_neg[row.names(div_classes_neg)%in% row.names(div_classes_pos),]
) 

# Sunburst plot
pdf(file="plots/ms2_bothpolarities_classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(bothPolarities_classes), bothPolarities_classes$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(bothPolarities_classes), "Number of spectra"=bothPolarities_classes$frequency), file="plots/ms2_bothpolarities_classes_sunburst.csv", row.names=FALSE)

# ---------- Sunburst plot of classes of all samples ----------
# Merge pos and neg div_classes
div_classes <- merge(div_classes_pos, div_classes_neg, by="row.names", all.x=TRUE, all.y=TRUE)
colnames(div_classes) <- c("superclass", "div_classes_pos", "div_classes_neg")
div_classes[is.na(div_classes)] <- 0
div_classes$frequency <- apply(X=div_classes[,c(2:ncol(div_classes))], MARGIN=1, FUN=function(x) { sum(x) })
rownames(div_classes) <- div_classes[,1]

# Plot
pdf(file="plots/ms2_classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(div_classes), div_classes$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_classes), "Number of spectra"=div_classes$frequency), file="plots/ms2_classes_sunburst.csv", row.names=FALSE)


# ---------- Variable Selection ----------
# BORUTA
#library(doParallel)
#cl <- makeCluster(16)
#registerDoParallel(cl)
sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=class_list, sel_factor=mzml_pheno_growth_stress_samples, sel_colors=mzml_pheno_colors$cols, sel_method="best", tune_length=10, pvalue=0.01, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms2_classes_select_boruta_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
f.heatmap.selected_features(feat_list = class_list, sel_feat=sel_boruta_rf, filename="plots/ms2_classes_select_boruta.pdf", main="BORUTA")
#stopCluster(cl)



# ############################## MS2 Superclasses statistics ##############################



# ---------- Histogram ----------
pdf(file="plots/ms2_superclasses_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(unlist(superclass_list)))
dev.off()



# ---------- Variation partitioning ----------
mzml_pheno_growth_stress_samples <- as.factor(c("Stress", "Stress", "Stress", "Stress", "Stress", "Stress", 
												"Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth",
												"Stress", "Stress", "Stress", "Stress", "Stress", "Stress",
												"Growth", "Growth", "Growth", "Growth", "Growth", "Growth",
												"Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress"))
model_varpart <- varpart(superclass_list, ~ mzml_pheno_samples, ~ mzml_pheno_growth_stress_samples)

mzml_pheno_methanol_samples <- as.factor(c("1", "1", "1", "10", "10", "10", 
																					 "1", "1", "1", "10", "10", "10",
																					 "0", "0", "0",
																					 "1", "1", "1", "10", "10", "10", 
																					 "1", "1", "1", "10", "10", "10",
																					 "0", "0", "0",
																					 "1", "1", "1", "10", "10", "10"))
model_varpart <- varpart(scale(feat_list), ~ mzml_pheno_samples, ~ mzml_pheno_methanol_samples)

# Plot results
pdf(file="plots/ms2_superclasses_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","MeOH"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- Diversity ----------
# Plot unique features
pdf(paste("plots/ms2_superclasses_richness_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$unique ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Number of unique features", xlab="treatment", ylab="number of unique features")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$unique, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot total features
pdf(paste("plots/ms2_superclasses_richness_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$features ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Number of features", xlab="treatment", ylab="number of features")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$features, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Shannon index
pdf(paste("plots/ms2_superclasses_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$shannon ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Pielou evenness for species
pdf(paste("plots/ms2_superclasses_diversity_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$pielou ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Pielou\'s evenness", xlab="treatment", ylab="Pielou diversity index (J)")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$pielou, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Sunburst plot
pdf(file="plots/ms2_superclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_superclasses), div_superclasses$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_superclasses), "Number of spectra"=div_superclasses$frequency), file="plots/ms2_superclasses_sunburst.csv", row.names=FALSE)

# superclass level
onlyPos_superclasses <- data.frame(class = row.names(onlyPos_classes), frequency = onlyPos_classes$frequency)
onlyPos_superclasses <-	separate(onlyPos_superclasses, col = class, into = c("Organic", "Superclass", "Class"), sep = ";")
onlyPos_superclasses <- onlyPos_superclasses %>% 
	group_by(Superclass) %>% 
	summarise(frequency = sum(frequency))
onlyPos_superclasses <- as.data.frame(onlyPos_superclasses)
row.names(onlyPos_superclasses) <- onlyPos_superclasses$Superclass

onlyNeg_superclasses <- data.frame(class = row.names(onlyNeg_classes), frequency = onlyNeg_classes$frequency)
onlyNeg_superclasses <-	separate(onlyNeg_superclasses, col = class, into = c("Organic", "Superclass", "Class"), sep = ";")
onlyNeg_superclasses <- onlyNeg_superclasses %>% 
	group_by(Superclass) %>% 
	summarise(frequency = sum(frequency))
onlyNeg_superclasses <- as.data.frame(onlyNeg_superclasses)
row.names(onlyNeg_superclasses) <- onlyNeg_superclasses$Superclass

bothPolarities_superclasses <- data.frame(class = row.names(bothPolarities_classes), frequency = bothPolarities_classes$frequency)
bothPolarities_superclasses <-	separate(bothPolarities_superclasses, col = class, into = c("Organic", "Superclass", "Class"), sep = ";")
bothPolarities_superclasses <- bothPolarities_superclasses %>% 
	group_by(Superclass) %>% 
	summarise(frequency = sum(frequency))
bothPolarities_superclasses <- as.data.frame(bothPolarities_superclasses)
row.names(bothPolarities_superclasses) <- bothPolarities_superclasses$Superclass

# Sunburst plot
pdf(file="plots/ms2_onlyPositive_superclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(onlyPos_superclasses), onlyPos_superclasses$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(onlyPos_superclasses), "Number of spectra"=onlyPos_superclasses$frequency), file="plots/ms2_onlyPositive_superclasses_sunburst.csv", row.names=FALSE)



# ---------- Variable Selection ----------
# BORUTA
#library(doParallel)
#cl <- makeCluster(16)
#registerDoParallel(cl)
sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=superclass_list, sel_factor=mzml_pheno_growth_stress_samples, sel_colors=mzml_pheno_colors$cols, sel_method="best", tune_length=10, pvalue=0.01, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms2_superclasses_select_boruta_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
f.heatmap.selected_features(feat_list=superclass_list, sel_feat=sel_boruta_rf, filename="plots/ms2_superclasses_select_boruta.pdf", main="BORUTA")
#stopCluster(cl)



#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################

library(pdftools)

pdf_convert("plots/ms1_diversity_shannon.pdf", 
						format = "png", 
						filenames = "plots/ms1_diversity_shannon.png", 
						dpi = 200)
