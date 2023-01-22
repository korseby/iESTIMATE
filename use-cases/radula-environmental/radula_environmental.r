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
#nSlaves <- 24
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
setwd("C:/Users/kaitl/Desktop/Research/UNB/Research Project/Germany/Bioinformatics/Radula_environmental")
#setwd("/Users/kristian/Desktop/Projekte/Habilitation/Mosses/radula_environmental")

# Data directory
mzml_dir <- paste(getwd(),"raw",sep="/")

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



# ---------- Load iESTIMATE functions ----------
source("_functions.r")



# ---------- Load meta-data ----------
radula_environmental_factors <- read.table(file=paste0("Sample-Information-environmental.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")


# ---------- Peak detection ----------
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

# comp_list
comp_list <- cbind(feat_list_pos[, c(rownames(ms1_def_pos)[ms1_def_pos$has_ms2==1])],
				   feat_list_neg[, c(rownames(ms1_def_neg)[ms1_def_neg$has_ms2==1])] )
colnames(comp_list) <- c(paste0(rownames(ms1_def_pos)[ms1_def_pos$has_ms2==1], "_pos"),
						 paste0(rownames(ms1_def_neg)[ms1_def_neg$has_ms2==1], "_neg") )

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
class_list <- cbind(class_list_neg)
rownames(class_list) <- gsub(x=rownames(class_list_neg), pattern="\\.auto.*", replacement="")
#class_list <- cbind(sapply(unique(colnames(class_list)[duplicated(colnames(class_list))]), function(x) rowSums(class_list[,grepl(paste(x, "$", sep=""), colnames(class_list))])), class_list[,!duplicated(colnames(class_list)) & ! duplicated(colnames(class_list), fromLast=TRUE)])

# superclass_list
superclass_list <- cbind(superclass_list_neg)
rownames(superclass_list) <- gsub(x=rownames(superclass_list), pattern="\\.auto.*", replacement="")
#superclass_list <- cbind(sapply(unique(colnames(superclass_list)[duplicated(colnames(superclass_list))]), function(x) rowSums(superclass_list[,grepl(paste(x, "$", sep=""), colnames(superclass_list))])), superclass_list[,!duplicated(colnames(superclass_list)) & ! duplicated(colnames(superclass_list), fromLast=TRUE)])

# div_class and div_superclass
library(dplyr)
div_classes <- rbind(div_classes_pos, div_classes_neg)

div_classes <- div_classes %>%
	group_by(rownames(div_classes)) %>%
	summarise(
		frequency = sum(frequency)
	)
div_classes <- as.data.frame(div_classes)

rm(classifiers_class)

# Search for MS1 spectra that have useable MS2 spectra
ms2_names_pos <- paste(names(ms2_spectra_pos), "_pos", sep = "")
ms2_names_neg <- paste(names(ms2_spectra_neg), "_neg", sep = "")
ms2_combined_names <- c(ms2_names_neg, ms2_names_pos)
mod_feat_list <- feat_list[, ms2_combined_names]

# make Canada specific stuff
mod_feat_list_can <- mod_feat_list[1:9, ]
month_factor_can <- as.factor(c("August", "August", "August", "June", "June", "June", "May", "May", "May"))
superclass_list_neg_can <- superclass_list_neg[1:9, ]
superclass_list_neg_can_flavonoids <- cbind(superclass_list_neg_can[, 73:77], superclass_list_neg_can[, 68], superclass_list_neg_can[, 80:82])

# Canada Flavonoid boxplot
can_flavonoids <- as.data.frame(superclass_list_neg_can_flavonoids)
can_flavonoids <- rowSums(can_flavonoids)
df_flavonoids <- data.frame(
	Abundance = can_flavonoids,
	Month = c("August", "August", "August", "June", "June", "June", "May", "May", "May")
) 
boxplot(df_flavonoids$Abundance ~ df_flavonoids$Month, 
				xlab = "Month", ylab = "Flavonoids Detected",
				col = c("blue3", "green3", "yellow3"))
summary(aov(df_flavonoids$Abundance ~ df_flavonoids$Month))
TukeyHSD(aov(df_flavonoids$Abundance ~ df_flavonoids$Month))


# ############################## MS1 statistics ##############################



# ---------- Histogram ----------
pdf(file="plots/ms1_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(feat_list))
dev.off()

# ---------- Variation partitioning ----------
# Plot results
pdf(file="plots/ms1_varpart2.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
model_varpart <- varpart(scale(feat_list), ~ radula_environmental_factors$East, ~radula_environmental_factors$North, ~radula_environmental_factors$Location, as.numeric(as.factor(radula_environmental_factors$Sample.name)))
plot(model_varpart, Xnames=c("Longtitude", "Latitude", "Country", "Samples"), cutoff=0, cex=1, id.size=1.2, digits=1, bg=c("blue","green", "yellow", "cyan"))

model_varpart <- varpart(scale(feat_list), ~ radula_environmental_factors$Location, ~radula_environmental_factors$Host.tree, ~radula_environmental_factors$Collection.Month, as.numeric(as.factor(radula_environmental_factors$Sample.name)))
plot(model_varpart, Xnames=c("Country", "Host Tree", "Month", "Samples"), cutoff=0, cex=1, id.size=1.2, digits=1, bg=c("yellow", "blue","green", "cyan"))

model_varpart <- varpart(scale(feat_list), ~ radula_environmental_factors$Location, ~ radula_environmental_factors$Host.tree.family, as.numeric(as.factor(radula_environmental_factors$Collection.Month)))
plot(model_varpart, Xnames=c("Location","Host tree family","Collection month"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green","red"))
dev.off()



# ---------- PCA ----------
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

pdf(file="plots/ms1_pca_23.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca <- prcomp(feat_list, center=TRUE)
plot(ms1_pca$x[, 2], ms1_pca$x[,3], pch=19, main="PCA: Grouping of samples",
		 xlab=paste0("PC2: ", format(summary(ms1_pca)$importance[2, 2] * 100, digits=3), " % variance"),
		 ylab=paste0("PC3: ", format(summary(ms1_pca)$importance[2, 3] * 100, digits=3), " % variance"),
		 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca$x[, 2], ms1_pca$x[,3], labels=names(feat_list[, 5]), col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()



# ==================== Host tree ====================
# ---------- Diversity ----------
mzml_pheno_samples <- radula_environmental_factors$Host.tree
mzml_pheno_colors <- data.frame(cols=2:(length(unique(mzml_pheno_samples))+1))
rownames(mzml_pheno_colors) <- unique(mzml_pheno_samples)
mzml_pheno_colors_samples <- sapply(mzml_pheno_samples, function(x) { x <- mzml_pheno_colors[,1][which(x==unique(mzml_pheno_samples))] } )

# Plot unique features
pdf(paste("plots/ms1_host_tree_richness_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$unique ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Number of unique features", xlab="treatment", ylab="number of unique features")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$unique, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot total features
pdf(paste("plots/ms1_host_tree_richness_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$features ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Number of features", xlab="treatment", ylab="number of features")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$features, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Shannon index
pdf(paste("plots/ms1_host_tree_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$shannon ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Pielou evenness for species
pdf(paste("plots/ms1_host_tree_diversity_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$pielou ~ mzml_pheno_samples, col=mzml_pheno_colors$cols, names=NA, main="Pielou\'s evenness", xlab="treatment", ylab="Pielou diversity index (J)")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$pielou, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# ---------- Environmental Variables: Compound Classes --------

radula_environmental_factors$Sample.name <- gsub("-", ".", radula_environmental_factors$Sample.name)
radula_environmental_factors$Host.tree <- gsub("(Ostrya virginiana and Betula papyrifera)", "Betulaceae", radula_environmental_factors$Host.tree)

#All superclasses
superclass_list_df <- as.data.frame(superclass_list)
superclass_list_df2 <- cbind(superclass_list_df, radula_environmental_factors)

heatmap.2(x=as.matrix(superclass_list), cexRow=0.6, cexCol=0.5,
					col=colorRampPalette(c("violet", "darkblue", "darkgreen", "yellow"), alpha=100, bias=1)(100), breaks = seq(0,10,0.1),
					trace="none", margins=c(max(nchar(colnames(class_list)))/5, max(nchar(rownames(class_list)))/3),
					key=TRUE, key.title="Color key", density.info='density', denscol="black")

#Specialized Metabolites
superclass_list_special <- cbind(superclass_list[, 186:187],superclass_list[134], 
																 superclass_list[61:62], superclass_list[58:59],
																 superclass_list[50], superclass_list[48], 
																 superclass_list[46], superclass_list[44], 
																 superclass_list[39], superclass_list[28], 
																 superclass_list[41])
superclass_list_special_df <- as.data.frame(superclass_list_special)

#Flavonoids
superclass_list_special_flavonoid <- cbind(superclass_list[, 149], 
																					 superclass_list[164:168],
																					 superclass_list[172:180],
																					 superclass_list[185])
superclass_list_special_flavonoid_df <- as.data.frame(superclass_list_special_flavonoid)

heatmap.2(x=as.matrix(superclass_list_special_df), cexRow=0.6, cexCol=0.8,
					#col=colorRampPalette(c(brewer.pal(n = 11, name = "RdBu")), alpha=1, bias=1)(256)
					col=colorRampPalette(c(brewer.pal(n = 11, name = "RdBu")), alpha=100, bias=1)(200), breaks = seq(0,20,0.1),
					trace="none", margins=c(max(nchar(colnames(superclass_list_special_df)))/3, max(nchar(rownames(superclass_list_special_df)))/2),
					key=TRUE, key.title="Color key", density.info='density', denscol="black")

#Flavonoid boxplot
superclass_list_flavonoids <- cbind(superclass_list[, 73:77], superclass_list[, 68])
flavonoids <- as.data.frame(superclass_list_flavonoids)
flavonoids <- rowSums(flavonoids)
flavonoids <- c(flavonoids[19:27], flavonoids[10:18], flavonoids[7:9], flavonoids[4:6], flavonoids[1:3])
df_flavonoids <- data.frame(
	Features = flavonoids,
	Sample = labels(flavonoids)
) 
df_flavonoids <- cbind(df_flavonoids, radula_environmental_factors)

boxplot(df_flavonoids$Features ~ df_flavonoids$Host.tree, 
				xlab = "Host Tree", ylab = "Flavonoids Detected",
				col = rainbow(7))
summary(aov(df_flavonoids$Features ~ df_flavonoids$Collection.Month))
TukeyHSD(aov(df_flavonoids$Features ~ df_flavonoids$Collection.Month))

#Flavonoid glycoside boxplot
superclass_list_flavonoids_glycoside <- cbind(superclass_list[, 80:82])
flavonoids_glycoside <- as.data.frame(superclass_list_flavonoids_glycoside)
flavonoids_glycoside <- rowSums(flavonoids_glycoside)
df_flavonoids_glycoside <- data.frame(
	Features = flavonoids_glycoside,
	Sample = labels(flavonoids_glycoside)
) 
df_flavonoids_glycoside <- cbind(df_flavonoids_glycoside, radula_environmental_factors)

boxplot(df_flavonoids_glycoside$Features ~ df_flavonoids_glycoside$Host.tree, 
				xlab = "Host Tree", ylab = "Flavonoid Glycosides Detected",
				col = rainbow(7))
summary(aov(df_flavonoids_glycoside$Features ~ df_flavonoids_glycoside$Collection.Month))
TukeyHSD(aov(df_flavonoids_glycoside$Features ~ df_flavonoids_glycoside$Collection.Month))

# Terpenoid boxplot
superclass_list_terpenoids <- cbind(superclass_list[, 29:30], superclass_list[, 32], superclass_list[, 34:35])
terpenoids <- as.data.frame(superclass_list_terpenoids)
terpenoids <- rowSums(terpenoids)
df_terpenoids <- data.frame(
	Features = terpenoids,
	Sample = labels(terpenoids)
) 
df_terpenoids <- cbind(df_terpenoids, radula_environmental_factors)

boxplot(df_terpenoids$Features ~ df_terpenoids$Host.tree, 
				xlab = "Host Tree", ylab = "Terpenoids Detected",
				col = rainbow(7))
summary(aov(df_terpenoids$Features ~ df_terpenoids$Collection.Month))
TukeyHSD(aov(df_terpenoids$Features ~ df_terpenoids$Collection.Month))

# Stilbene boxplot
superclass_list_stilbenes <- cbind(superclass_list[, 84])
stilbenes <- as.data.frame(superclass_list_stilbenes)
stilbenes <- rowSums(stilbenes)
df_stilbenes <- data.frame(
	Features = stilbenes,
	Sample = labels(stilbenes)
) 
df_stilbenes <- cbind(df_stilbenes, radula_environmental_factors)

boxplot(df_stilbenes$Features ~ df_stilbenes$Host.tree, 
				xlab = "Host Tree", ylab = "Stilbenes Detected",
				col = rainbow(7))
summary(aov(df_stilbenes$Features ~ df_stilbenes$Collection.Month))
TukeyHSD(aov(df_stilbenes$Features ~ df_stilbenes$Collection.Month))

# Sugar boxplot
superclass_list_sugars <- cbind(superclass_list[, 53], superclass_list[, 55], superclass_list[, 59:61])
sugars <- as.data.frame(superclass_list_sugars)
sugars <- rowSums(sugars)
df_sugars <- data.frame(
	Features = sugars,
	Sample = labels(sugars)
) 
df_sugars <- cbind(df_sugars, radula_environmental_factors)

boxplot(df_sugars$Features ~ df_sugars$Location, 
				xlab = "Location", ylab = "Sugars Detected",
				col = rainbow(7))
summary(aov(df_sugars$Features ~ df_sugars$Collection.Month))
TukeyHSD(aov(df_sugars$Features ~ df_sugars$Collection.Month))

# Fatty acids
superclass_list_fatty_acids <- cbind(superclass_list[, 25:28])
fatty_acids <- as.data.frame(superclass_list_fatty_acids)
fatty_acids <- rowSums(fatty_acids)
df_fatty_acids <- data.frame(
	Features = fatty_acids,
	Sample = labels(fatty_acids)
) 
df_fatty_acids <- cbind(df_fatty_acids, radula_environmental_factors)

boxplot(df_fatty_acids$Features ~ df_fatty_acids$Location, 
				xlab = "Location", ylab = "Fatty Acids Detected",
				col = rainbow(7))
summary(aov(df_fatty_acids$Features ~ df_fatty_acids$Collection.Month))
TukeyHSD(aov(df_fatty_acids$Features ~ df_fatty_acids$Collection.Month))

# Fancy boxplot
library(ggpubr)
library(viridis)
library(ggplot2)
b <- ggplot(class_list_df, aes(y = Features, x = Collection.Month, fill = Location)) + scale_fill_viridis(discrete = TRUE, option = "D")
b2 <- b + geom_hline(yintercept=0) + geom_boxplot() + theme_bw() 
b2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 20))
#my_comparisons <- list( c("AA1", "AA10"), c("JA1", "JA10"), c("SA1", "SA10"), c("NAA1", "NAA10"), c("BAP1", "BAP10"))
#b3 <- b2 + stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(0, 0, 0, 0, 0)) 
#b3

# ---------- Variable Selection ----------
# Random Forest
sel_rf <- f.select_features_random_forest(feat_matrix=feat_list, sel_factor=as.factor(mzml_pheno_samples), sel_colors=mzml_pheno_colors$cols, tune_length=10, quantile_threshold=0.99, plot_roc_filename="plots/ms1_host_tree_select_rf_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_rf$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=feat_list, sel_feat=sel_rf$`_selected_variables_`, filename="plots/ms1_host_tree_select_rf.pdf", main="Random Forest")

# BORUTA
#library(doParallel)
#cl <- makeCluster(24)
#registerDoParallel(cl)
#sel_boruta_rf <- f.select_features_random_forest_boruta(feat_matrix=feat_list, sel_factor=as.factor(mzml_pheno_samples), sel_colors=mzml_pheno_colors$cols, sel_method="best", tune_length=2, pvalue=0.01, mcAdj=TRUE, maxRuns=10000, plot_roc_filename="plots/ms1_select_host_tree_boruta_roc.pdf")
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_boruta_rf[["_selected_variables_"]])))
#f.heatmap.selected_features(feat_list=feat_list, sel_feat=sel_boruta_rf, filename="plots/ms1_host_tree_select_boruta.pdf", main="BORUTA")
#stopCluster(cl)

#make factors
host_tree_factors <- as.factor(radula_environmental_factors$Host.tree)
location_factors <- as.factor(radula_environmental_factors$Location)
month_factors <- as.factor(radula_environmental_factors$Collection.Month)
host_tree_family_factors <- as.factor(radula_environmental_factors$Host.tree.family)

# Missing value imputation, NAs set as 0 above
feat_list[feat_list==0]<-0.1
mod_feat_list[mod_feat_list==0]<-0.1

#[which(is.na(feat_list))] <- median(na.omit(as.numeric(unlist(feat_list))))
#feat_list[which(is.na(feat_list))] <- runif(length(which(is.na(feat_list))), min=0, max=0.0000001)
#feat_list <- scale(feat_list, scale=FALSE, center=FALSE)
#feat_list[is.na(feat_list)] <- 0
#feat_list[which(feat_list < 0)] <- 0
#feat_list[is.infinite(feat_list)] <- 0
#feat_list <- feat_list[!apply(feat_list, MARGIN=1, function(x) max(x,na.rm=TRUE) == min(x,na.rm=TRUE)),]

# PLS (doesn't like 0s in data)
sel_pls <- f.select_features_pls(feat_matrix=mod_feat_list, sel_factor=location_factors, sel_colors=mzml_pheno_colors$cols, components=2, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/ms1_mod_select_pls_roc_location4.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=mod_feat_list, sel_feat=sel_pls$`_selected_variables_`, filename="plots/ms1_mod_select_pls_location4.pdf", main="PLS")

# components = 2 is fine for location, host tree might have more (6)

# PLS (class level)
sel_pls <- f.select_features_pls(feat_matrix=class_list, sel_factor=host_tree_family_factors, sel_colors=mzml_pheno_colors$cols, components=4, tune_length=10, quantile_threshold=0.99, plot_roc_filename="plots/ms1_class_mod_select_pls_roc_host_tree_family.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=class_list, sel_feat=sel_pls$`_selected_variables_`, filename="plots/ms1_class_mod_select_pls_host_tree_family.pdf", main="PLS")

# PLS (superclass level)
sel_pls <- f.select_features_pls(feat_matrix=superclass_list, sel_factor=host_tree_factors, sel_colors=mzml_pheno_colors$cols, components=4, tune_length=10, quantile_threshold=0.99, plot_roc_filename="plots/ms1_superclass_mod_select_pls_roc_tree.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=superclass_list, sel_feat=sel_pls$`_selected_variables_`, filename="plots/ms1_superclass_mod_select_pls_tree.pdf", main="PLS")

# sPLS-DA
# First: Variation partitioning
model_varpart <- varpart(scale(mod_feat_list), ~ location_factors, ~ location_factors)

# 10% of features correlate with factor
model_varpart_corr <- trunc(model_varpart$part$indfract$Adj.R.squared[2] * ncol(mod_feat_list))
model_splsda_keepx <- trunc(seq(model_varpart_corr, 1,length.out=length(unique(location_factors))))

sel_splsda <- f.select_features_splsda(feat_matrix=mod_feat_list, sel_colors=mzml_pheno_colors$cols, sel_factor=location_factors, tune_components=length(unique(location_factors)) - 1, sel_components=2, folds_number=8, keepx=model_splsda_keepx, plot_roc_filename="plots/ms1_mod_select_splsda_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_splsda$'_selected_variables_')))
f.heatmap.selected_features(feat_list=mod_feat_list, sel_feat=sel_splsda$'_selected_variables_', filename="plots/ms1_mod_select_splsda.pdf", main="sPLS-DA")

# re-load functions (for tweaking)
source("_functions.r")

# convert PDFs to PNG
library(pdftools)

pdf_convert("plots/ms1_mod_select_pls_location4.pdf", 
						format = "png", 
						page = 1,
						filenames = "plots/mod_select_pls_location4.png", 
						dpi = 200)

#feat_list, sel_feat, sel_names=NULL, sample_colors=NULL, filename, main, scale="col", plot_width=6, plot_height=5, cex_col=0.5, cex_row=0

# prepare MetaboLights MAF
#Location
results_location_table <- read.table(file="data/results_location.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)
results_location_xcms <- gsub(x=results_location_table$Features, pattern="_.*", replacement="")
results_location_mode <- gsub(x=results_location_table$Features, pattern=".*_", replacement="")
results_location_smiles <- results_location_table$SMILES
results_location_name <- results_location_table$Identification
for (i in 1:length(results_location_smiles)) {
	if ((results_location_smiles[i] == "-") | (results_location_smiles[i] == "n.d") | (results_location_smiles[i] == "")) {
		results_location_smiles[i] <- NA
	}
	if ((results_location_name[i] == "-") | (results_location_name[i] == "n.d") | (results_location_name[i] == "")) {
		results_location_name[i] <- NA
	}
	if (is.na(results_location_smiles[i])) {
		results_location_xcms[i] <- NA
		results_location_mode[i] <- NA
	}
}

results_location_xcms <- na.omit(results_location_xcms)
results_location_mode <- na.omit(results_location_mode)
results_location_smiles <- na.omit(results_location_smiles)
results_location_name <- na.omit(results_location_name)

#host_tree_species
results_host_tree_species_table <- read.table(file="data/results_host_tree_species.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)
results_host_tree_species_xcms <- gsub(x=results_host_tree_species_table$Features, pattern="_.*", replacement="")
results_host_tree_species_mode <- gsub(x=results_host_tree_species_table$Features, pattern=".*_", replacement="")
results_host_tree_species_smiles <- results_host_tree_species_table$SMILES
results_host_tree_species_name <- results_host_tree_species_table$Identification
for (i in 1:length(results_host_tree_species_smiles)) {
	if ((results_host_tree_species_smiles[i] == "-") | (results_host_tree_species_smiles[i] == "n.d") | (results_host_tree_species_smiles[i] == "")) {
		results_host_tree_species_smiles[i] <- NA
	}
	if ((results_host_tree_species_name[i] == "-") | (results_host_tree_species_name[i] == "n.d") | (results_host_tree_species_name[i] == "")) {
		results_host_tree_species_name[i] <- NA
	}
	if (is.na(results_host_tree_species_smiles[i])) {
		results_host_tree_species_xcms[i] <- NA
		results_host_tree_species_mode[i] <- NA
	}
}

results_host_tree_species_xcms <- na.omit(results_host_tree_species_xcms)
results_host_tree_species_mode <- na.omit(results_host_tree_species_mode)
results_host_tree_species_smiles <- na.omit(results_host_tree_species_smiles)
results_host_tree_species_name <- na.omit(results_host_tree_species_name)

#host_tree_family
results_host_tree_family_table <- read.table(file="data/results_host_tree_family.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)
results_host_tree_family_xcms <- gsub(x=results_host_tree_family_table$Features, pattern="_.*", replacement="")
results_host_tree_family_mode <- gsub(x=results_host_tree_family_table$Features, pattern=".*_", replacement="")
results_host_tree_family_smiles <- results_host_tree_family_table$SMILES
results_host_tree_family_name <- results_host_tree_family_table$Identification
for (i in 1:length(results_host_tree_family_smiles)) {
	if ((results_host_tree_family_smiles[i] == "-") | (results_host_tree_family_smiles[i] == "n.d") | (results_host_tree_family_smiles[i] == "")) {
		results_host_tree_family_smiles[i] <- NA
	}
	if ((results_host_tree_family_name[i] == "-") | (results_host_tree_family_name[i] == "n.d") | (results_host_tree_family_name[i] == "")) {
		results_host_tree_family_name[i] <- NA
	}
	if (is.na(results_host_tree_family_smiles[i])) {
		results_host_tree_family_xcms[i] <- NA
		results_host_tree_family_mode[i] <- NA
	}
}

results_host_tree_family_xcms <- na.omit(results_host_tree_family_xcms)
results_host_tree_family_mode <- na.omit(results_host_tree_family_mode)
results_host_tree_family_smiles <- na.omit(results_host_tree_family_smiles)
results_host_tree_family_name <- na.omit(results_host_tree_family_name)

#Combine identifications
results_xcms <- c(results_location_xcms, results_host_tree_species_xcms, results_host_tree_family_xcms)
results_mode <- c(results_location_mode, results_host_tree_species_mode, results_host_tree_family_mode)
results_smiles <- c(results_location_smiles, results_host_tree_species_smiles, results_host_tree_family_smiles)
results_name <- c(results_location_name, results_host_tree_species_name, results_host_tree_family_name)

ms1_def_neg$smiles <- ""
ms1_def_neg[which(rownames(ms1_def_neg) %in% results_xcms), "smiles"] <- results_smiles

ms1_def_pos$smiles <- ""

# Annotate feat_list in neg-mode
results_maf_neg <- ms1_def_neg
f.export_maf(results_maf_neg, "data/results_maf_neg.tsv")
f.annotate_maf_classes(maf_input="data/results_maf_neg.tsv", maf_output="data/results_maf_neg_annotated_classes.tsv")
f.annotate_maf_compounds(maf_input="data/results_maf_neg_annotated_classes.tsv", maf_output="data/results_maf_neg_annotated_compounds.tsv", polarity="neg", xcms_id=results_xcms, pol_mode=results_mode, smiles=results_smiles, names=results_name)

# Annotate feat_list in pos-mode
results_maf_pos <- ms1_def_pos
f.export_maf(results_maf_pos, "data/results_maf_pos.tsv")
f.annotate_maf_classes(maf_input="data/results_maf_pos.tsv", maf_output="data/results_maf_pos_annotated_classes.tsv")
#f.annotate_maf_compounds(maf_input="data/ingroup_maf_pos_annotated_classes.tsv", maf_output="data/ingroup_maf_pos_annotated_compounds.tsv", polarity="pos", xcms_id=riccia_ingroup_xcms, pol_mode=riccia_ingroup_mode, smiles=riccia_ingroup_smiles, names=riccia_ingroup_name)

