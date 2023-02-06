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
setwd("~/Desktop/Projekte/Habilitation/Mosses/riccia_chemotax/")

# Data directory
mzml_dir <- "~/Desktop/Projekte/Habilitation/Mosses/riccia_chemotax/"

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
SIRIUS_VERSION <- 5                  # Version of SIRIUS
mzabs <- 0.01                        # Absolute mass error (in seconds) used for merging MS/MS spectra
mzppm <- 35                          # ppm error used for merging MS/MS spectra
msms.intensity.threshold <- 10       # Minimum intensity value for MS/MS peaks

# Preparations for plotting
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)

# Save and load RData
#save.image()
load(".RData")



# ---------- Peak detection ----------
source("_functions.r")

setMSnbaseFastLoad(TRUE)

source("peak_detection_pos.r")
source("peak_detection_neg.r")



# ---------- Create merged pos+neg objects ----------
# factors
mzml_names <- gsub(x=mzml_names_pos, pattern="\\.pos.*", replacement="")
mzml_pheno_samples <- mzml_pheno_samples_pos
mzml_pheno_colors <- mzml_pheno_colors_pos
mzml_pheno_colors_samples <- mzml_pheno_colors_samples_pos

# feat_list
feat_list <- cbind(feat_list_pos, feat_list_neg)
rownames(feat_list) <- gsub(x=rownames(feat_list_pos), pattern="\\.auto.*", replacement="")
colnames(feat_list) <- c(paste0(colnames(feat_list_pos),"_pos"), paste0(colnames(feat_list_neg),"_neg"))

# comp_list
comp_list <- cbind(feat_list_pos[, c(rownames(ms1_def_pos)[which(ms1_def_pos$has_ms2==1)])],
				   feat_list_neg[, c(rownames(ms1_def_neg)[which(ms1_def_neg$has_ms2==1)])] )
colnames(comp_list) <- c(paste0(rownames(ms1_def_pos)[which(ms1_def_pos$has_ms2==1)], "_pos"),
						 paste0(rownames(ms1_def_neg)[which(ms1_def_neg$has_ms2==1)], "_neg") )

# bina_list
bina_list <- cbind(bina_list_pos[, c(rownames(ms1_def_pos)[which(ms1_def_pos$has_ms2==1)])],
				   bina_list_neg[, c(rownames(ms1_def_neg)[which(ms1_def_neg$has_ms2==1)])] )
rownames(bina_list) <- gsub(x=rownames(bina_list_pos), pattern="\\.auto.*", replacement="")
colnames(bina_list) <- c(paste0(rownames(ms1_def_pos)[which(ms1_def_pos$has_ms2==1)], "_pos"),
						 paste0(rownames(ms1_def_neg)[which(ms1_def_neg$has_ms2==1)], "_neg") )

# uniq_list
uniq_list <- cbind(uniq_list_pos[, c(rownames(ms1_def_pos)[which(ms1_def_pos$has_ms2==1)])],
				   uniq_list_neg[, c(rownames(ms1_def_neg)[which(ms1_def_neg$has_ms2==1)])] )
rownames(uniq_list) <- gsub(x=rownames(uniq_list_pos), pattern="\\.auto.*", replacement="")
colnames(uniq_list) <- c(paste0(rownames(ms1_def_pos)[which(ms1_def_pos$has_ms2==1)], "_pos"),
						 paste0(rownames(ms1_def_neg)[which(ms1_def_neg$has_ms2==1)], "_neg") )

# Create data frame
model_div             <- data.frame(features=apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div$richness    <- apply(X=bina_list, MARGIN=1, FUN=function(x) { sum(x) } )
#model_div$menhinick   <- apply(X=bina_list, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div$shannon     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div$pielou      <- apply(X=scale(comp_list, center=F, scale=T), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list, species), index="chao")
model_div$simpson     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div$inverse     <- apply(X=comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div$fisher      <- apply(X=comp_list, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div$unique      <- apply(X=uniq_list, MARGIN=1, FUN=function(x) { sum(x) })

# Remove NAs if present
model_div[is.na(model_div)] <- 0

# class_list
class_list <- cbind(class_list_pos, class_list_neg)
rownames(class_list) <- gsub(x=rownames(class_list_pos), pattern="\\.auto.*", replacement="")
class_list <- cbind(sapply(unique(colnames(class_list)[duplicated(colnames(class_list))]), function(x) rowSums(class_list[,grepl(paste(x, "$", sep=""), colnames(class_list))])), class_list[,!duplicated(colnames(class_list)) & ! duplicated(colnames(class_list), fromLast=TRUE)])

# class_int_list
class_int_list <- cbind(class_int_list_pos, class_int_list_neg)
rownames(class_int_list) <- gsub(x=rownames(class_int_list_pos), pattern="\\.auto.*", replacement="")
class_int_list <- cbind(sapply(unique(colnames(class_int_list)[duplicated(colnames(class_int_list))]), function(x) rowSums(class_int_list[,grepl(paste(x, "$", sep=""), colnames(class_int_list))])), class_int_list[,!duplicated(colnames(class_int_list)) & ! duplicated(colnames(class_int_list), fromLast=TRUE)])

# superclass_list
superclass_list <- cbind(superclass_list_pos, superclass_list_neg)
rownames(superclass_list) <- gsub(x=rownames(superclass_list_pos), pattern="\\.auto.*", replacement="")
superclass_list <- cbind(sapply(unique(colnames(superclass_list)[duplicated(colnames(superclass_list))]), function(x) rowSums(superclass_list[,grepl(paste(x, "$", sep=""), colnames(superclass_list))])), superclass_list[,!duplicated(colnames(superclass_list)) & ! duplicated(colnames(superclass_list), fromLast=TRUE)])

# superclass_int_list
superclass_int_list <- cbind(superclass_int_list_pos, superclass_int_list_neg)
rownames(superclass_int_list) <- gsub(x=rownames(superclass_int_list_pos), pattern="\\.auto.*", replacement="")
superclass_int_list <- cbind(sapply(unique(colnames(superclass_int_list)[duplicated(colnames(superclass_int_list))]), function(x) rowSums(superclass_int_list[,grepl(paste(x, "$", sep=""), colnames(superclass_int_list))])), superclass_int_list[,!duplicated(colnames(superclass_int_list)) & ! duplicated(colnames(superclass_int_list), fromLast=TRUE)])



# ############################## MS1 statistics ##############################



# ---------- Histogram ----------
pdf(file="plots/ms1_comp_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(comp_list))
dev.off()



# ---------- Variation partitioning ----------
model_varpart <- varpart(scale(comp_list), ~ mzml_pheno_samples, ~ mzml_pheno_samples)

# Plot results
pdf(file="plots/ms1_comp_list_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca <- prcomp(comp_list, scale=TRUE, center=TRUE)

# PC-Axes 1+2
pdf(paste("plots/ms1_comp_list_pca12.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 1], model_pca$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC1: ", format(summary(model_pca)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples, bg=mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 1], model_pca$x[,2], labels=mzml_names, col=mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()

# PC-Axes 2+3
pdf(paste("plots/ms1_comp_list_pca23.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 2], model_pca$x[,3], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples, bg=mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 2], model_pca$x[,3], labels=mzml_names, col=mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()

# 3D PCA
library(rgl)
rgl.open() # Open a new RGL device
bg3d("white")
#par3d(windowRect = c(100, 100, 1200, 1200))
pca_3dplot <- plot3d(model_pca$x[, 1:3], col=mzml_pheno_colors_samples, theta=30, phi=15)
text3d(model_pca$x[, 1:3], texts=mzml_names, col=mzml_pheno_colors_samples)
rgl.postscript("___.pdf", "pdf")

htmlwidgets::saveWidget(file="____.html", rglwidget(width=1600,height=1600), libdir="libs",	selfcontained=TRUE)

#rgl.bg(color = "grey87")
#pca4<-NULL
#pca4<-pca3$scores[1:nrow(data_client),]

#pca3d(as.matrix(pca4),components = 1:3,group=data_client[1:nrow(data_client),]$cluster,show.ellipses=TRUE,show.axes=TRUE,show.axe.titles = TRUE,axes.color = "#404040",title="PCA",
#	  show.shapes=TRUE,show.shadows = FALSE,show.plane=FALSE,fancy=FALSE,show.group.labels = FALSE,show.centroids = FALSE,
#	  labels.col=data_client[1:nrow(data_client),]$RColor,ellipse.ci=0.95,col=data_client[1:nrow(data_client),]$RColor,radius=1.4,shape="sphere")
#rgl.viewpoint( theta = -58, phi = 20, fov = 60, zoom = 0.55, 
#			   scale = par3d("scale"), interactive = TRUE, 
#			   type = c("userviewpoint", "modelviewpoint"))

#Saving plot
#snapshotPCA3d(file="snapshotPCA3d.bmp")
#library(pca3d)
#gr <- factor(comp_list)
#pca3d(model_pca, group=gr)



# ---------- Diversity ----------
# Plot unique features
pdf(paste("plots/ms1_comp_list_richness_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$unique ~ mzml_pheno_samples, col=mzml_pheno_colors, names=NA, main="Number of unique compounds", xlab="treatment", ylab="number of unique compounds")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$unique, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot total features
pdf(paste("plots/ms1_comp_list_richness_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$features ~ mzml_pheno_samples, col=mzml_pheno_colors, names=NA, main="Number of compounds", xlab="treatment", ylab="number of compounds")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$features, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Shannon index
pdf(paste("plots/ms1_comp_list_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$shannon ~ mzml_pheno_samples, col=mzml_pheno_colors, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$shannon, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Pielou evenness for species
pdf(paste("plots/ms1_comp_list_diversity_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div$pielou ~ mzml_pheno_samples, col=mzml_pheno_colors, names=NA, main="Pielou\'s evenness", xlab="treatment", ylab="Pielou diversity index (J)")
text(1:length(unique(mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div$pielou, term=as.factor(mzml_pheno_samples))
text(1:length(unique(mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()



# ---------- Variable Selection ----------
# Load identification
riccia_id <- read.table(file="data/riccia_id.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)
riccia_id$sel_names <- riccia_id$XCMS.name

for (i in 1:length(riccia_id$sel_names)) {
	if ((riccia_id$Compound.Name[i] != "-") & (riccia_id$Compound.Name[i] != "n.d")) {
		riccia_id$sel_names[i] <- substr(riccia_id$Compound.Name[i], 1, 50)
	}
}

# PLS
sel_pls <- f.select_features_pls(feat_matrix=comp_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=2, tune_length=10, quantile_threshold=0.975, plot_roc_filename="plots/ms1_comp_list_select_pls_roc.pdf")
print(paste("Number of selected compounds:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=comp_list, sel_feat=sel_pls$`_selected_variables_`, sel_names=riccia_id$sel_names, plot_width=8, plot_height=4, sample_colors=mzml_pheno_colors_samples, filename="plots/ms1_comp_list_select_pls.pdf", main="")
#f.heatmap.selected_features(feat_list=comp_list, sel_feat=sel_pls$`_selected_variables_`, sel_names=paste0("     ",sel_pls$`_selected_variables_`), plot_width=8, plot_height=3.3, sample_colors=mzml_pheno_colors_samples, filename="plots/ms1_comp_list_select_pls.pdf", main="")
#f.heatmap.selected_features(feat_list=comp_list, sel_feat=sel_pls$`_selected_variables_`, sel_names=paste0("     ",sel_pls$`_selected_variables_`), sample_colors=mzml_pheno_colors_samples, filename="plots/ms1_comp_list_select_pls.pdf", main="PLS")
sel_pls$`_selected_variables_`
sel_pls$`_model_r2_`
sel_pls$`_multiclass_metrics_`
#boxplot(comp_list[,"FT0024_pos"] ~ mzml_pheno_samples, col=mzml_pheno_colors)

# PLS with bina_list
sel_pls <- f.select_features_pls(feat_matrix=bina_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=2, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/ms1_bina_list_select_pls_roc.pdf")
print(paste("Number of selected compounds:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=bina_list, sel_feat=sel_pls$`_selected_variables_`, sel_names=paste0("     ",sel_pls$`_selected_variables_`), scale="none", sample_colors=mzml_pheno_colors_samples, filename="plots/ms1_bina_list_select_pls.pdf", main="PLS")



# ############################## MS2 classes statistics ##############################



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



# ---------- Histogram ----------
pdf(file="plots/ms2_classes_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(unlist(class_list)), main="Histogram", xlab="class_list")
dev.off()



# ---------- Variation partitioning ----------
model_varpart <- varpart(scale(class_list), ~ mzml_pheno_samples, ~ mzml_pheno_samples)

# Plot results
pdf(file="plots/ms2_class_list_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca <- prcomp(class_list, scale=TRUE, center=TRUE)

# PC-Axes 1+2
pdf(paste("plots/ms2_class_list_pca12.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 1], model_pca$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC1: ", format(summary(model_pca)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples, bg=mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 1], model_pca$x[,2], labels=mzml_names, col=mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()

# PC-Axes 2+3
pdf(paste("plots/ms2_class_list_pca23.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 2], model_pca$x[,3], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples, bg=mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 2], model_pca$x[,3], labels=mzml_names, col=mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()



# ---------- Variable Selection ----------
# PLS
sel_pls <- f.select_features_pls(feat_matrix=class_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=2, tune_length=10, quantile_threshold=0.975, plot_roc_filename="plots/ms2_class_list_select_pls_roc.pdf")
print(paste("Number of selected compounds:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=class_list, sel_feat=sel_pls$`_selected_variables_`, sample_colors=mzml_pheno_colors_samples, plot_width=4, plot_height=4, filename="plots/ms2_class_list_select_pls.pdf", main="")
sel_pls$`_multiclass_metrics_`
sel_pls$`_model_r2_`



# ---------- Histogram ----------
pdf(file="plots/ms2_classes_int_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(unlist(class_int_list)), main="Histogram", xlab="class_int_list")
dev.off()



# ---------- Variation partitioning ----------
model_varpart <- varpart(class_int_list, ~ mzml_pheno_samples, ~ mzml_pheno_samples)

# Plot results
pdf(file="plots/ms2_class_int_list_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca <- prcomp(class_int_list, scale=FALSE, center=TRUE)

# PC-Axes 1+2
pdf(paste("plots/ms2_class_int_list_pca12.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 1], model_pca$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC1: ", format(summary(model_pca)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples, bg=mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 1], model_pca$x[,2], labels=mzml_names, col=mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()

# PC-Axes 2+3
pdf(paste("plots/ms2_class_int_list_pca23.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 2], model_pca$x[,3], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples, bg=mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 2], model_pca$x[,3], labels=mzml_names, col=mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()



# ---------- Variable Selection ----------
# PLS
sel_pls <- f.select_features_pls(feat_matrix=class_int_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=2, tune_length=10, quantile_threshold=0.975, plot_roc_filename="plots/ms2_class_int_list_select_pls_roc.pdf")
print(paste("Number of selected compounds:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=class_int_list, sel_feat=sel_pls$`_selected_variables_`, sample_colors=mzml_pheno_colors_samples, filename="plots/ms2_class_int_list_select_pls.pdf", main="PLS")
sel_pls$`_multiclass_metrics_`
sel_pls$`_model_r2_`



# ############################## MS2 superclasses statistics ##############################



# ---------- Sunburst plot of superclasses of all samples ----------
# Merge pos and neg div_superclasses
div_superclasses <- merge(div_superclasses_pos, div_superclasses_neg, by="row.names", all.x=TRUE, all.y=TRUE)
colnames(div_superclasses) <- c("superclass", "div_superclasses_pos", "div_superclasses_neg")
div_superclasses[is.na(div_superclasses)] <- 0
div_superclasses$frequency <- apply(X=div_superclasses[,c(2:ncol(div_superclasses))], MARGIN=1, FUN=function(x) { sum(x) })
rownames(div_superclasses) <- div_superclasses[,1]

# Plot
pdf(file="plots/ms2_superclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(div_superclasses), div_superclasses$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound superclasses"=rownames(div_superclasses), "Number of spectra"=div_superclasses$frequency), file="plots/ms2_superclasses_sunburst.csv", row.names=FALSE)



# ---------- Histogram ----------
pdf(file="plots/ms2_superclasses_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(unlist(superclass_list)), main="Histogram", xlab="superclass_list")
dev.off()



# ---------- Variation partitioning ----------
model_varpart <- varpart(scale(superclass_list), ~ mzml_pheno_samples, ~ mzml_pheno_samples)

# Plot results
pdf(file="plots/ms2_superclass_list_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca <- prcomp(superclass_list, scale=TRUE, center=TRUE)

# PC-Axes 1+2
pdf(paste("plots/ms2_superclass_list_pca12.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 1], model_pca$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC1: ", format(summary(model_pca)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples, bg=mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 1], model_pca$x[,2], labels=mzml_names, col=mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()

# PC-Axes 2+3
pdf(paste("plots/ms2_superclass_list_pca23.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 2], model_pca$x[,3], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples, bg=mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 2], model_pca$x[,3], labels=mzml_names, col=mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()



# ---------- Variable Selection ----------
# PLS
sel_pls <- f.select_features_pls(feat_matrix=superclass_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=2, tune_length=10, quantile_threshold=0.925, plot_roc_filename="plots/ms2_superclass_list_select_pls_roc.pdf")
print(paste("Number of selected compounds:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=superclass_list, sel_feat=sel_pls$`_selected_variables_`, sample_colors=mzml_pheno_colors_samples, plot_width=4, plot_height=4, filename="plots/ms2_superclass_list_select_pls.pdf", main="")
sel_pls$`_multiclass_metrics_`
sel_pls$`_model_r2_`



# ---------- Histogram ----------
pdf(file="plots/ms2_superclasses_int_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(unlist(superclass_int_list)), main="Histogram", xlab="superclass_int_list")
dev.off()



# ---------- Variation partitioning ----------
model_varpart <- varpart(scale(superclass_int_list), ~ mzml_pheno_samples, ~ mzml_pheno_samples)

# Plot results
pdf(file="plots/ms2_superclass_int_list_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca <- prcomp(superclass_int_list, scale=TRUE, center=TRUE)

# PC-Axes 1+2
pdf(paste("plots/ms2_superclass_int_list_pca12.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 1], model_pca$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC1: ", format(summary(model_pca)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples, bg=mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 1], model_pca$x[,2], labels=mzml_names, col=mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()

# PC-Axes 2+3
pdf(paste("plots/ms2_superclass_int_list_pca23.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 2], model_pca$x[,3], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples, bg=mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 2], model_pca$x[,3], labels=mzml_names, col=mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()



# ---------- Variable Selection ----------
# PLS
sel_pls <- f.select_features_pls(feat_matrix=superclass_int_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=3, tune_length=10, quantile_threshold=0.925, plot_roc_filename="plots/ms2_superclass_int_list_select_pls_roc.pdf")
print(paste("Number of selected compounds:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=superclass_int_list, sel_feat=sel_pls$`_selected_variables_`, sample_colors=mzml_pheno_colors_samples, filename="plots/ms2_superclass_int_list_select_pls.pdf", main="PLS")
sel_pls$`_multiclass_metrics_`
sel_pls$`_model_r2_`



# ############################## Molecular Descriptors ##############################



# ---------- Read SIRIUS identification ----------
# Get SMILES
sirius_compound_identifications_pos <- read.table(file=paste0("data/pos_ms2_compound_identifications.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
sirius_compound_identifications_neg <- read.table(file=paste0("data/neg_ms2_compound_identifications.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
sirius_compounds <- c(sirius_compound_identifications_pos$smiles, sirius_compound_identifications_neg$smiles)
sirius_compound_ids <- c(gsub(x=sirius_compound_identifications_pos$id, pattern=".*_", replacement=""), gsub(x=sirius_compound_identifications_neg$id, pattern=".*_", replacement=""))

ms1_def_pos$smiles <- ""
for (i in 1:nrow(sirius_compound_identifications_pos)) {
	ms1_def_pos$smiles[which(rownames(ms1_def_pos)==gsub(x=sirius_compound_identifications_pos$id[i], pattern=".*_", replacement=""))] <- sirius_compound_identifications_pos$smiles[i]
}

ms1_def_neg$smiles <- ""
for (i in 1:nrow(sirius_compound_identifications_neg)) {
	ms1_def_neg$smiles[which(rownames(ms1_def_neg)==gsub(x=sirius_compound_identifications_neg$id[i], pattern=".*_", replacement=""))] <- sirius_compound_identifications_neg$smiles[i]
}



# ---------- Calculate molecular descriptors with RDKit ----------
# Write list of SMILES for rdkit
writeLines(sirius_compounds, "data/compound_identifications_smiles.txt")

# Execute rdkit script
system2(command="./rdkit_calculate_descriptors.py", args=c("-i","data/compound_identifications_smiles.txt", "-o","data/descriptors_rdkit.csv"), wait=TRUE)

# Read rdkit_calculate_descriptors.py output
rdkit_descriptors <- read.csv(file="data/descriptors_rdkit.csv", header=TRUE, sep=";", fill=TRUE, dec=".", stringsAsFactors=FALSE, check.names=FALSE)
rownames(rdkit_descriptors) <- colnames(sirius_compounds)
rdkit_descriptors[,2:ncol(rdkit_descriptors)] <- as.data.frame(apply(X=rdkit_descriptors[,2:ncol(rdkit_descriptors)], MARGIN=2, FUN=function(x) { as.numeric(x) }))
rdkit_descriptors <- rdkit_descriptors[,2:ncol(rdkit_descriptors)]

# Properly format NAs
for (i in 1:nrow(rdkit_descriptors)) {
	rdkit_descriptors[i, ] <- as.numeric(rdkit_descriptors[i, ])
}
for (i in 1:nrow(rdkit_descriptors)) {
	rdkit_descriptors[i, as.character(rdkit_descriptors[i,]) %in% 'list(NULL)'] <- 0
	rdkit_descriptors[i, as.character(rdkit_descriptors[i,]) %in% 'NA'] <- 0
	rdkit_descriptors[i, as.character(rdkit_descriptors[i,]) %in% 'NULL'] <- 0
	rdkit_descriptors[i, as.character(rdkit_descriptors[i,]) %in% 'NaN'] <- 0
	rdkit_descriptors[i, as.character(rdkit_descriptors[i,]) %in% 'Inf'] <- 0
}
rdkit_descriptors <- as.data.frame(rdkit_descriptors)



# ---------- Calculate molecular descriptors from COCONUT library ----------
# Query descriptors for MetFrag results in positive mode
library(httr)
library(jsonlite)

coconut_descriptors <- NULL
for (smile in sirius_compounds) {
	if ( ! is.na(nchar(smile)) ) {
		coconut_descriptors[[smile]] <- plyr::rbind.fill(coconut_descriptors[[smile]], tryCatch(f.query_coconut(smile), error=function(x) data.frame(smiles=smile)))
		if (length(coconut_descriptors[[smile]]) < 1) coconut_descriptors[[smile]] <- NULL
	}
}

# Calculate general molecular descriptors from SMILES
number_of_phosphorus <- unlist(lapply(X=lapply(X=coconut_descriptors, FUN=function(x) { x <- x$molecular_formula }), FUN=function(x) { x <- as.numeric(unlist(Map(function(x) { x <- gsub("(?<=[A-z])(?![0-9])","1",x,perl = TRUE); table(factor(rep(gsub("\\d+", "", x),as.numeric(gsub("\\D+", "", x))), levels="P"))}, regmatches(x, gregexpr("[A-z]+(\\d+)?", x))))) } ))

c_n_ratio <- unlist(lapply(X=coconut_descriptors, FUN=function(x) { x <- x$number_of_carbons })) / unlist(lapply(X=coconut_descriptors, FUN=function(x) { x <- x$number_of_nitrogens }))
c_n_ratio[is.infinite(c_n_ratio)] <- 0

c_n_p_ratio <- c_n_ratio / number_of_phosphorus
c_n_p_ratio[is.infinite(c_n_p_ratio)] <- 0
c_n_p_ratio[is.na(c_n_p_ratio)] <- 0

c_p_ratio <- unlist(lapply(X=coconut_descriptors, FUN=function(x) { x <- x$number_of_carbons })) / number_of_phosphorus
c_p_ratio[is.infinite(c_p_ratio)] <- 0

for (i in 1:length(coconut_descriptors)) {
	coconut_descriptors[[i]]$number_of_phosphorus <- number_of_phosphorus[i]
	coconut_descriptors[[i]]$c_n_ratio <- c_n_ratio[i]
	coconut_descriptors[[i]]$c_n_p_ratio <- c_n_p_ratio[i]
	coconut_descriptors[[i]]$c_p_ratio <- c_p_ratio[i]
}

# Build data frame
descriptors <- NULL
for (i in 1:length(coconut_descriptors)) {
	descriptors <- plyr::rbind.fill(descriptors, as.data.frame(coconut_descriptors[[i]][1, ! grepl(x=colnames(coconut_descriptors[[i]]), pattern="(fragments|ertl|absolute_smiles)", perl=TRUE)]))
}

# Filter data frame with only the relevant columns
coconut_properties <- c("npl_noh_score", "npl_score", "npl_sugar_score", "number_of_carbons",
						"number_of_nitrogens", "number_of_oxygens", "number_of_rings", "max_number_of_rings", "min_number_of_rings",
						"total_atom_number", "bond_count",
						"alogp", "alogp2", "amralogp", "apol", "bpol", "eccentricConnectivityIndexDescriptor",
						"fmfDescriptor", "fsp3", "fragmentComplexityDescriptor", "hBondAcceptorCount", "hBondDonorCount",
						"hybridizationRatioDescriptor", "kappaShapeIndex1", "kappaShapeIndex2", "kappaShapeIndex3",
						"manholdlogp", "petitjeanNumber", "petitjeanShapeTopo", "petitjeanShapeGeom", "lipinskiRuleOf5Failures",
						"numberSpiroAtoms", "vabcDescriptor", "vertexAdjMagnitude",
						"zagrebIndex", "tpsaEfficiency", "weinerPolarityNumber",
						"pfCounts.count",
						"number_of_phosphorus", "c_n_ratio", "c_n_p_ratio", "c_p_ratio")
descriptors <- descriptors[, coconut_properties]

# Properly format NAs
for (i in 1:nrow(descriptors)) {
	descriptors[i, as.character(descriptors[i,]) %in% 'list(NULL)'] <- 0
	descriptors[i, as.character(descriptors[i,]) %in% 'NA'] <- 0
	descriptors[i, as.character(descriptors[i,]) %in% 'NULL'] <- 0
	descriptors[i, as.character(descriptors[i,]) %in% 'NaN'] <- 0
}

# Remove rows with all NA or all 0
#descriptors <- descriptors[rowSums(is.na(descriptors)) != ncol(descriptors), ]
#descriptors <- na.omit(descriptors)
#descriptors <- descriptors[(rowSums(descriptors) != 0), ]

# Convert to numeric
for (i in 1:nrow(descriptors)) {
	descriptors[i, ] <- as.numeric(unlist(descriptors[i,]))
}



# ---------- Calculate molecular descriptors with CDK ----------
library(rcdk)

cdk_descriptors <- NULL
for (i in sirius_compounds) {
	# Get Structure from SMILES
	cdk_mol = parse.smiles(i)[[1]]
	
	# Get simple measures
	cdk_atoms = get.atoms(cdk_mol)
	cdk_bonds = get.bonds(cdk_mol)
	
	# Calculate simple measures
	cdk_num_atoms = as.factor(unlist(lapply(cdk_atoms, get.symbol)))
	cdk_num_atoms = tapply(cdk_num_atoms, cdk_num_atoms, length)
	numC = as.numeric(cdk_num_atoms["C"])
	numN = as.numeric(cdk_num_atoms["N"])
	numP = as.numeric(cdk_num_atoms["P"])
	numO = as.numeric(cdk_num_atoms["O"])
	CNRatio = as.numeric(numC / numN)
	
	# Calculate descriptors and restrict to only "constitutional" and "topological"
	cdk_mol_des_cats = get.desc.categories()
	cdk_mol_des_names = c(get.desc.names(cdk_mol_des_cats[3]), get.desc.names(cdk_mol_des_cats[4]))
	cdk_mol_des = as.data.frame(eval.desc(cdk_mol, cdk_mol_des_names))
	cdk_descriptors <- plyr::rbind.fill(cdk_descriptors, cbind(data.frame(numC=numC, numN=numN, numP=numP, numO=numO, CNRatio=CNRatio), cdk_mol_des))
}

# Properly format NAs and convert to numeric
for (i in 1:nrow(cdk_descriptors)) {
	cdk_descriptors[i, as.character(cdk_descriptors[i,]) %in% 'list(NULL)'] <- 0
	cdk_descriptors[i, as.character(cdk_descriptors[i,]) %in% 'NA'] <- 0
	cdk_descriptors[i, as.character(cdk_descriptors[i,]) %in% 'NULL'] <- 0
	cdk_descriptors[i, as.character(cdk_descriptors[i,]) %in% 'NaN'] <- 0
	cdk_descriptors[i, ] <- as.numeric(unlist(cdk_descriptors[i,]))
}



# ---------- Build descriptor list ----------
# Table L: samples x metabolites
mdes_tab_l <- class_list#comp_list

# Table R: samples x species
mdes_tab_r <- as.data.frame.matrix(table(rownames(mdes_tab_l), mzml_pheno_samples))
rownames(mdes_tab_r) <- rownames(mdes_tab_l)

# Table Q: metabolites x traits
mdes_tab_q <- cdk_descriptors
mdes_tab_q[is.na(mdes_tab_q)] <- 0
#mdes_tab_l <- comp_list[, which(colnames(comp_list) %in% c(paste0(rownames(ms1_def_pos)[which(ms1_def_pos$smiles != "")], "_pos"), paste0(rownames(ms1_def_neg)[which(ms1_def_neg$smiles != "")], "_neg")))]
mdes_tab_l <- class_list[, which(colnames(class_list) %in% gsub(x=ms1_def_pos$primary_class,pattern=".*; ",replacement="")[which(gsub(x=ms1_def_pos$primary_class,pattern=".*; ",replacement="") != "")]), which(colnames(class_list) %in% gsub(x=ms1_def_neg$primary_class,pattern=".*; ",replacement="")[which(gsub(x=ms1_def_neg$primary_class,pattern=".*; ",replacement="") != "")]) ]
mdes_tab_l <- as.data.frame(mdes_tab_l)

# Perform matrix operation
mdes_list <- as.data.frame(as.matrix(mdes_tab_l) %*% as.matrix(mdes_tab_q))



# ---------- Perform Full dbRDA ----------
attach(mdes_list)
mdes_rda_formula <- formula(paste0("~ 0 + ", paste0(colnames(mdes_list),collapse=" + ")))
mdes_rda_y <- data.frame(model.matrix(mdes_rda_formula))
detach(mdes_list)

# Calculate overlay of descriptors on comp_list distances
model_mdes_dbrda <- vegan::dbrda(formula=comp_list ~ ., data=mdes_rda_y, distance="euclidean")
model_mdes_dbrda_scores <- vegan::scores(model_mdes_dbrda, choices=c(1,2))
model_mdes_dbrda_ef_descriptors <- envfit(model_mdes_dbrda, mdes_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
model_mdes_dbrda_fit_descriptors <- data.frame(r2=c(model_mdes_dbrda_ef_descriptors$vectors$r,model_mdes_dbrda_ef_descriptors$factors$r),
											   pvals=c(model_mdes_dbrda_ef_descriptors$vectors$pvals,model_mdes_dbrda_ef_descriptors$factors$pvals) )
rownames(model_mdes_dbrda_fit_descriptors) <- c(names(model_mdes_dbrda_ef_descriptors$vectors$r),names(model_mdes_dbrda_ef_descriptors$factors$r))
model_mdes_dbrda_fit_descriptors
write.csv(model_mdes_dbrda_fit_descriptors, file="plots/descriptors_dbrda_fit.csv", row.names=TRUE)


# Re-calculate dbRDA with selected variables
mdes_sel_list <- as.data.frame(as.matrix(mdes_tab_l) %*% as.matrix(mdes_tab_q[, rownames(model_mdes_dbrda_fit_descriptors)[which(model_mdes_dbrda_fit_descriptors$pvals<0.01)]]))
attach(mdes_sel_list)
mdes_rda_sel_formula <- formula(paste0("~ 0 + ", paste0(colnames(mdes_sel_list),collapse=" + ")))
mdes_rda_sel_y <- data.frame(model.matrix(mdes_rda_sel_formula))
detach(mdes_sel_list)
model_mdes_dbrda_sel <- vegan::dbrda(formula=comp_list ~ ., data=mdes_rda_sel_y, distance="euclidean")
model_mdes_dbrda_sel_scores <- vegan::scores(model_mdes_dbrda_sel, choices=c(1,2))

# Calculate overlay of descriptors on metabolite profile distances
model_mdes_dbrda_sel_ef_descriptors <- envfit(model_mdes_dbrda_sel, mdes_sel_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
model_mdes_dbrda_sel_fit_descriptors <- data.frame(r2=c(model_mdes_dbrda_sel_ef_descriptors$vectors$r,model_mdes_dbrda_sel_ef_descriptors$factors$r),
											   pvals=c(model_mdes_dbrda_sel_ef_descriptors$vectors$pvals,model_mdes_dbrda_sel_ef_descriptors$factors$pvals) )
rownames(model_mdes_dbrda_sel_fit_descriptors) <- c(names(model_mdes_dbrda_sel_ef_descriptors$vectors$r),names(model_mdes_dbrda_sel_ef_descriptors$factors$r))
model_mdes_dbrda_sel_fit_descriptors
write.csv(model_mdes_dbrda_sel_fit_descriptors, file="plots/descriptors_dbrda_sel_fit.csv", row.names=TRUE)


# Plot dbRDA for descriptors
pdf(file="plots/descriptors_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_mdes_dbrda_sel_scores$sites[,1])-1, max(model_mdes_dbrda_sel_scores$sites[,1])+1),
	 ylim=c(min(model_mdes_dbrda_sel_scores$sites[,2]), max(model_mdes_dbrda_sel_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_mdes_dbrda_sel)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_mdes_dbrda_sel)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_mdes_dbrda_sel, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_mdes_dbrda_sel, display="sites", labels=rownames(comp_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_mdes_dbrda_sel_ef_descriptors, cex=0.75, p.max=1, col="black")
dev.off()


# Calculate overlay of comp_list on metabolite profile distances
model_mdes_dbrda_ef_comp_list <- envfit(model_mdes_dbrda, comp_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
model_mdes_dbrda_fit_comp_list <- data.frame(r2=c(model_mdes_dbrda_ef_comp_list$vectors$r,model_mdes_dbrda_ef_comp_list$factors$r),
											 pvals=c(model_mdes_dbrda_ef_comp_list$vectors$pvals,model_mdes_dbrda_ef_comp_list$factors$pvals) )
rownames(model_mdes_dbrda_fit_comp_list) <- c(names(model_mdes_dbrda_ef_comp_list$vectors$r),names(model_mdes_dbrda_ef_comp_list$factors$r))
model_mdes_dbrda_fit_comp_list
write.csv(model_mdes_dbrda_fit_comp_list, file="plots/comp_list_dbrda_fit.csv", row.names=TRUE)

# Plot dbRDA for descriptors
pdf(file="plots/comp_list_dbrda_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_mdes_dbrda_scores$sites[,1])-1, max(model_mdes_dbrda_scores$sites[,1])+1),
	 ylim=c(min(model_mdes_dbrda_scores$sites[,2]), max(model_mdes_dbrda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_mdes_dbrda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_mdes_dbrda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_mdes_dbrda, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_mdes_dbrda, display="sites", labels=rownames(comp_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_mdes_dbrda_ef_comp_list, cex=0.75, p.max=1, col="black")
dev.off()


# Calculate overlay of class_list on metabolite profile distances
model_mdes_dbrda_ef_class_list <- envfit(model_mdes_dbrda, class_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
model_mdes_dbrda_fit_class_list <- data.frame(r2=c(model_mdes_dbrda_ef_class_list$vectors$r,model_mdes_dbrda_ef_class_list$factors$r),
											  pvals=c(model_mdes_dbrda_ef_class_list$vectors$pvals,model_mdes_dbrda_ef_class_list$factors$pvals) )
rownames(model_mdes_dbrda_fit_class_list) <- c(names(model_mdes_dbrda_ef_class_list$vectors$r),names(model_mdes_dbrda_ef_class_list$factors$r))
model_mdes_dbrda_fit_class_list
write.csv(model_mdes_dbrda_fit_class_list, file="plots/class_list_dbrda_fit.csv", row.names=TRUE)

# Plot dbRDA for descriptors
pdf(file="plots/class_list_dbrda_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_mdes_dbrda_scores$sites[,1])-1, max(model_mdes_dbrda_scores$sites[,1])+1),
	 ylim=c(min(model_mdes_dbrda_scores$sites[,2]), max(model_mdes_dbrda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_mdes_dbrda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_mdes_dbrda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_mdes_dbrda, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_mdes_dbrda, display="sites", labels=rownames(comp_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_mdes_dbrda_ef_class_list, cex=0.75, p.max=1, col="black")
dev.off()


# Calculate overlay of superclass_list on metabolite profile distances
model_mdes_dbrda_ef_superclass_list <- envfit(model_mdes_dbrda, superclass_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
model_mdes_dbrda_fit_superclass_list <- data.frame(r2=c(model_mdes_dbrda_ef_superclass_list$vectors$r,model_mdes_dbrda_ef_superclass_list$factors$r),
												   pvals=c(model_mdes_dbrda_ef_superclass_list$vectors$pvals,model_mdes_dbrda_ef_superclass_list$factors$pvals) )
rownames(model_mdes_dbrda_fit_superclass_list) <- c(names(model_mdes_dbrda_ef_superclass_list$vectors$r),names(model_mdes_dbrda_ef_superclass_list$factors$r))
model_mdes_dbrda_fit_superclass_list
write.csv(model_mdes_dbrda_fit_superclass_list, file="plots/superclass_list_dbrda_fit.csv", row.names=TRUE)

# Plot dbRDA for descriptors
pdf(file="plots/superclass_list_dbrda_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_mdes_dbrda_scores$sites[,1])-1, max(model_mdes_dbrda_scores$sites[,1])+1),
	 ylim=c(min(model_mdes_dbrda_scores$sites[,2]), max(model_mdes_dbrda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_mdes_dbrda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_mdes_dbrda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_mdes_dbrda, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_mdes_dbrda, display="sites", labels=rownames(comp_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_mdes_dbrda_ef_superclass_list, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- dbRDA of mdes_list on comp_list with stepwise forward selection ----------
attach(mdes_list)
mdes_rda_comp_formula <- formula(paste0("~ 0 + ", paste0(colnames(mdes_list),collapse=" + ")))
mdes_rda_comp_y <- data.frame(model.matrix(mdes_rda_comp_formula))
detach(mdes_list)

# Model with intercept only
model_0_mdes_rda_comp <- vegan::dbrda(formula=comp_list ~ 1, data=mdes_rda_comp_y, distance="euclidean")

# Model with all explanatory variables
model_1_mdes_rda_comp <- vegan::dbrda(formula=comp_list ~ ., data=mdes_rda_comp_y, distance="euclidean")

# Stepwise forward selection
model_step_mdes_rda_comp <- ordistep(object=model_0_mdes_rda_comp, scope=formula(model_1_mdes_rda_comp), Pin=0.1, Pout=0.5, direction="forward", perm.max=100000)
model_step_mdes_rda_comp_scores <- vegan::scores(model_step_mdes_rda_comp)

# dbRDA with selected model by permutation tests in constrained ordination
model_mdes_rda_comp <- vegan::dbrda(formula=as.formula(model_step_mdes_rda_comp$terms), data=mdes_rda_comp_y, distance="euclidean")
model_mdes_rda_comp_ef_formula <- update(as.formula(model_step_mdes_rda_comp$terms), model_mdes_rda_comp ~ .)
model_mdes_rda_comp_ef_factors <- as.factor(sapply(strsplit(as.character(model_mdes_rda_comp_ef_formula)[[3]], "\\+"), function(x) { x <- gsub("(\\`|^ | $)","",x) }))
model_mdes_rda_comp_ef <- envfit(formula=model_mdes_rda_comp_ef_formula, data=mdes_rda_comp_y, perm=10000)
model_mdes_rda_comp_scores <- vegan::scores(model_mdes_rda_comp)

# Goodness of fit statistic: Squared correlation coefficient
model_mdes_rda_comp_fit <- data.frame(r2=c(model_mdes_rda_comp_ef$vectors$r,model_mdes_rda_comp_ef$factors$r),
									  pvals=c(model_mdes_rda_comp_ef$vectors$pvals,model_mdes_rda_comp_ef$factors$pvals) )
rownames(model_mdes_rda_comp_fit) <- c(names(model_mdes_rda_comp_ef$vectors$r),names(model_mdes_rda_comp_ef$factors$r))
model_mdes_rda_comp_fit
write.csv(model_mdes_rda_comp_fit, file="plots/ms1_comp_list_dbrda_sel_fit.csv", row.names=TRUE)

# Plot results
pdf(file="plots/ms1_comp_list_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_mdes_rda_comp_scores$sites[,1])-1, max(model_mdes_rda_comp_scores$sites[,1])+1),
	 ylim=c(min(model_mdes_rda_comp_scores$sites[,2]), max(model_mdes_rda_comp_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_mdes_rda_comp)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_mdes_rda_comp)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_mdes_rda_comp, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_mdes_rda_comp, display="sites", labels=rownames(comp_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_mdes_rda_comp_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- dbRDA of mdes_list on class_list with stepwise forward selection ----------
attach(mdes_list)
mdes_rda_class_formula <- formula(paste0("~ 0 + ", paste0(colnames(mdes_list),collapse=" + ")))
mdes_rda_class_y <- data.frame(model.matrix(mdes_rda_class_formula))
detach(mdes_list)

# Model with intercept only
model_0_mdes_rda_class <- vegan::dbrda(formula=class_list ~ 1, data=mdes_rda_class_y, distance="euclidean")

# Model with all explanatory variables
model_1_mdes_rda_class <- vegan::dbrda(formula=class_list ~ ., data=mdes_rda_class_y, distance="euclidean")

# Stepwise forward selection
model_step_mdes_rda_class <- ordistep(object=model_0_mdes_rda_class, scope=formula(model_1_mdes_rda_class), Pin=0.1, Pout=0.5, direction="forward", perm.max=100000)
model_step_mdes_rda_class_scores <- vegan::scores(model_step_mdes_rda_class)

# dbRDA with selected model by permutation tests in constrained ordination
model_mdes_rda_class <- vegan::dbrda(formula=as.formula(model_step_mdes_rda_class$terms), data=mdes_rda_class_y, distance="euclidean")
model_mdes_rda_class_ef_formula <- update(as.formula(model_step_mdes_rda_class$terms), model_mdes_rda_class ~ .)
model_mdes_rda_class_ef_factors <- as.factor(sapply(strsplit(as.character(model_mdes_rda_class_ef_formula)[[3]], "\\+"), function(x) { x <- gsub("(\\`|^ | $)","",x) }))
model_mdes_rda_class_ef <- envfit(formula=model_mdes_rda_class_ef_formula, data=mdes_rda_class_y, perm=10000)
model_mdes_rda_class_scores <- vegan::scores(model_mdes_rda_class)

# Goodness of fit statistic: Squared correlation coefficient
model_mdes_rda_class_fit <- data.frame(r2=c(model_mdes_rda_class_ef$vectors$r,model_mdes_rda_class_ef$factors$r),
									   pvals=c(model_mdes_rda_class_ef$vectors$pvals,model_mdes_rda_class_ef$factors$pvals) )
rownames(model_mdes_rda_class_fit) <- c(names(model_mdes_rda_class_ef$vectors$r),names(model_mdes_rda_class_ef$factors$r))
model_mdes_rda_class_fit
write.csv(model_mdes_rda_class_fit, file="plots/ms2_class_list_dbrda_sel_fit.csv", row.names=TRUE)

# Plot results
pdf(file="plots/ms2_class_list_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_mdes_rda_class_scores$sites[,1])-1, max(model_mdes_rda_class_scores$sites[,1])+1),
	 ylim=c(min(model_mdes_rda_class_scores$sites[,2]), max(model_mdes_rda_class_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_mdes_rda_class)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_mdes_rda_class)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_mdes_rda_class, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_mdes_rda_class, display="sites", labels=rownames(class_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_mdes_rda_class_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- dbRDA of mdes_list on superclass_list with stepwise forward selection ----------
attach(mdes_list)
mdes_rda_superclass_formula <- formula(paste0("~ 0 + ", paste0(colnames(mdes_list),collapse=" + ")))
mdes_rda_superclass_y <- data.frame(model.matrix(mdes_rda_superclass_formula))
detach(mdes_list)

# Model with intercept only
model_0_mdes_rda_superclass <- vegan::dbrda(formula=superclass_list ~ 1, data=mdes_rda_superclass_y, distance="euclidean")

# Model with all explanatory variables
model_1_mdes_rda_superclass <- vegan::dbrda(formula=superclass_list ~ ., data=mdes_rda_superclass_y, distance="euclidean")

# Stepwise forward selection
model_step_mdes_rda_superclass <- ordistep(object=model_0_mdes_rda_superclass, scope=formula(model_1_mdes_rda_superclass), Pin=0.1, Pout=0.5, direction="forward", perm.max=100000)
model_step_mdes_rda_superclass_scores <- vegan::scores(model_step_mdes_rda_superclass)

# dbRDA with selected model by permutation tests in constrained ordination
model_mdes_rda_superclass <- vegan::dbrda(formula=as.formula(model_step_mdes_rda_superclass$terms), data=mdes_rda_superclass_y, distance="euclidean")
model_mdes_rda_superclass_ef_formula <- update(as.formula(model_step_mdes_rda_superclass$terms), model_mdes_rda_superclass ~ .)
model_mdes_rda_superclass_ef_factors <- as.factor(sapply(strsplit(as.character(model_mdes_rda_superclass_ef_formula)[[3]], "\\+"), function(x) { x <- gsub("(\\`|^ | $)","",x) }))
model_mdes_rda_superclass_ef <- envfit(formula=model_mdes_rda_superclass_ef_formula, data=mdes_rda_superclass_y, perm=10000)
model_mdes_rda_superclass_scores <- vegan::scores(model_mdes_rda_superclass)

# Goodness of fit statistic: Squared correlation coefficient
model_mdes_rda_superclass_fit <- data.frame(r2=c(model_mdes_rda_superclass_ef$vectors$r,model_mdes_rda_superclass_ef$factors$r),
											pvals=c(model_mdes_rda_superclass_ef$vectors$pvals,model_mdes_rda_superclass_ef$factors$pvals) )
rownames(model_mdes_rda_superclass_fit) <- c(names(model_mdes_rda_superclass_ef$vectors$r),names(model_mdes_rda_superclass_ef$factors$r))
model_mdes_rda_superclass_fit
write.csv(model_mdes_rda_superclass_fit, file="plots/ms2_superclass_list_dbrda_sel_fit.csv", row.names=TRUE)

# Plot results
pdf(file="plots/ms2_superclass_list_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_mdes_rda_superclass_scores$sites[,1])-1, max(model_mdes_rda_superclass_scores$sites[,1])+1),
	 ylim=c(min(model_mdes_rda_superclass_scores$sites[,2]), max(model_mdes_rda_superclass_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_mdes_rda_superclass)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_mdes_rda_superclass)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_mdes_rda_superclass, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_mdes_rda_superclass, display="sites", labels=rownames(superclass_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_mdes_rda_superclass_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- PLS ----------
# PLS
sel_pls <- f.select_features_pls(feat_matrix=mdes_list, sel_factor=mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=2, tune_length=10, quantile_threshold=0.995, plot_roc_filename="plots/descriptors_select_pls_roc.pdf")
print(paste("Number of selected descriptors:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=mdes_list, sel_feat=sel_pls$`_selected_variables_`, sample_colors=mzml_pheno_colors_samples, plot_width=4, plot_height=4, filename="plots/descriptors_select_pls.pdf", main="")
sel_pls$`_multiclass_metrics_`
sel_pls$`_model_r2_`



# ############################## Functional trait-image analysis ##############################



# ---------- Analyze habitus bioimages ----------
#BiocManager::install("EBImage")
library(EBImage)

# Read image files
img_files <- list.files("data/", pattern="*.jpg", recursive=FALSE, full.names=TRUE)
img_names <- gsub(x=img_files, pattern=".*/", replacement="")
#img_names <- gsub(x=img_names, pattern="_hab.*", replacement="")
riccia_img <- list()

# Get histogram and save values
for (i in 1:length(img_files)) {
	# Extraxt histogram
	riccia_img[[img_names[i]]] <- list()
	riccia_img[[img_names[i]]]$image <- readImage(img_files[i])
	riccia_img[[img_names[i]]]$hist <- hist(riccia_img[[img_names[i]]]$image)
	
	# Red
	riccia_img[[img_names[i]]]$red$x <- riccia_img[[img_names[i]]]$hist$red$mids[2:(length(riccia_img[[img_names[i]]]$hist$red$mids)-1)]
	riccia_img[[img_names[i]]]$red$y <- riccia_img[[img_names[i]]]$hist$red$counts[2:(length(riccia_img[[img_names[i]]]$hist$red$counts)-1)]
	riccia_img[[img_names[i]]]$red$loess_x <- loess.smooth(x=riccia_img[[img_names[i]]]$red$x, y=riccia_img[[img_names[i]]]$red$y, family="gaussian", span=0.1, col="red", lpars=list(col="red", lwd=2))$x
	riccia_img[[img_names[i]]]$red$loess_y <- loess.smooth(x=riccia_img[[img_names[i]]]$red$x, y=riccia_img[[img_names[i]]]$red$y, family="gaussian", span=0.1, col="red", lpars=list(col="red", lwd=2))$y
	riccia_img[[img_names[i]]]$red$max <- riccia_img[[img_names[i]]]$red$x[riccia_img[[img_names[i]]]$red$y==max(riccia_img[[img_names[i]]]$red$y)][1]
	
	# Green
	riccia_img[[img_names[i]]]$green$x <- riccia_img[[img_names[i]]]$hist$green$mids[2:(length(riccia_img[[img_names[i]]]$hist$green$mids)-1)]
	riccia_img[[img_names[i]]]$green$y <- riccia_img[[img_names[i]]]$hist$green$counts[2:(length(riccia_img[[img_names[i]]]$hist$green$counts)-1)]
	riccia_img[[img_names[i]]]$green$loess_x <- loess.smooth(x=riccia_img[[img_names[i]]]$green$x, y=riccia_img[[img_names[i]]]$green$y, family="gaussian", span=0.1, col="green", lpars=list(col="green", lwd=2))$x
	riccia_img[[img_names[i]]]$green$loess_y <- loess.smooth(x=riccia_img[[img_names[i]]]$green$x, y=riccia_img[[img_names[i]]]$green$y, family="gaussian", span=0.1, col="green", lpars=list(col="green", lwd=2))$y
	riccia_img[[img_names[i]]]$green$max <- riccia_img[[img_names[i]]]$green$x[riccia_img[[img_names[i]]]$green$y==max(riccia_img[[img_names[i]]]$green$y)][1]
	
	# Blue
	riccia_img[[img_names[i]]]$blue$x <- riccia_img[[img_names[i]]]$hist$blue$mids[2:(length(riccia_img[[img_names[i]]]$hist$blue$mids)-1)]
	riccia_img[[img_names[i]]]$blue$y <- riccia_img[[img_names[i]]]$hist$blue$counts[2:(length(riccia_img[[img_names[i]]]$hist$blue$counts)-1)]
	riccia_img[[img_names[i]]]$blue$loess_x <- loess.smooth(x=riccia_img[[img_names[i]]]$blue$x, y=riccia_img[[img_names[i]]]$blue$y, family="gaussian", span=0.1, col="blue", lpars=list(col="blue", lwd=2))$x
	riccia_img[[img_names[i]]]$blue$loess_y <- loess.smooth(x=riccia_img[[img_names[i]]]$blue$x, y=riccia_img[[img_names[i]]]$blue$y, family="gaussian", span=0.1, col="blue", lpars=list(col="blue", lwd=2))$y
	riccia_img[[img_names[i]]]$blue$max <- riccia_img[[img_names[i]]]$blue$x[riccia_img[[img_names[i]]]$blue$y==max(riccia_img[[img_names[i]]]$blue$y)][1]
}

# Plot histograms
pdf(file="plots/img_hist.pdf", encoding="ISOLatin1", pointsize=12, width=6, height=4, family="Helvetica")
for (i in 1:length(img_files)) {
	# Red
	plot(x=riccia_img[[img_names[i]]]$red$x, y=riccia_img[[img_names[i]]]$red$y, ylim=c(0,max(c(riccia_img[[img_names[i]]]$red$y, riccia_img[[img_names[i]]]$green$y), riccia_img[[img_names[i]]]$blue$y)), type="points", col="red", pch=19, xlab="Intensity", ylab="counts", main="")
	lines(x=riccia_img[[img_names[i]]]$red$loess_x, y=riccia_img[[img_names[i]]]$red$loess_y, col="red", lwd=2)
	
	# Green
	points(x=riccia_img[[img_names[i]]]$green$x, y=riccia_img[[img_names[i]]]$green$y, type="points", col="green", pch=19)
	lines(x=riccia_img[[img_names[i]]]$green$loess_x, y=riccia_img[[img_names[i]]]$green$loess_y, col="green", lwd=2)
	
	# Blue
	points(x=riccia_img[[img_names[i]]]$blue$x, y=riccia_img[[img_names[i]]]$blue$y, type="points", col="blue", pch=19)
	lines(x=riccia_img[[img_names[i]]]$blue$loess_x, y=riccia_img[[img_names[i]]]$blue$loess_y, col="blue", lwd=2)
	
	if (i==1) mtext(text="(a)", at=c(-0.12), line=2.5, font=2, cex=1.4)
	else if (i==2) mtext(text="(b)", at=c(-0.12), line=2.5, font=2, cex=1.4)
	else if (i==3) mtext(text="(c)", at=c(-0.12), line=2.5, font=2, cex=1.4)
}
dev.off()

# Extract max RGB values and save in data frame
img_trait_list <- data.frame(species=rep(img_names, each=3))
rownames(img_trait_list) <- mzml_names
for (i in 1:length(img_files)) {
	#img_trait_list[img_trait_list$species==img_names[i], "red"] <- riccia_img[[img_names[i]]]$red$max
	#img_trait_list[img_trait_list$species==img_names[i], "green"] <- riccia_img[[img_names[i]]]$green$max
	#img_trait_list[img_trait_list$species==img_names[i], "blue"] <- riccia_img[[img_names[i]]]$blue$max
	for (j in 1:length(riccia_img[[img_names[i]]]$red$y)) {
		img_trait_list[img_trait_list$species==img_names[i], paste0("red_",j)] <- riccia_img[[img_names[i]]]$red$y[j]
		img_trait_list[img_trait_list$species==img_names[i], paste0("green_",j)] <- riccia_img[[img_names[i]]]$green$y[j]
		img_trait_list[img_trait_list$species==img_names[i], paste0("bluered_",j)] <- riccia_img[[img_names[i]]]$blue$y[j]
	}
}
img_trait_list <- img_trait_list[, -c(1)]



# ---------- Full dbRDA of RGB x comp_list ----------
attach(as.data.frame(comp_list))
img_trait_rda_formula <- formula(paste0("~ 0 + ", paste0(colnames(comp_list),collapse=" + ")))
img_trait_rda_y <- data.frame(model.matrix(img_trait_rda_formula))
detach(as.data.frame(comp_list))

# Full dbRDA with compound list
model_img_trait_dbrda <- vegan::dbrda(formula=img_trait_list ~ ., data=img_trait_rda_y, distance="bray")
model_img_trait_dbrda_scores <- vegan::scores(model_img_trait_dbrda, choices=c(1,2))
model_img_trait_dbrda_ef <- envfit(model_img_trait_dbrda, comp_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
model_img_trait_dbrda_fit <- data.frame(r2=c(model_img_trait_dbrda_ef$vectors$r,model_img_trait_dbrda_ef$factors$r),
										pvals=c(model_img_trait_dbrda_ef$vectors$pvals,model_img_trait_dbrda_ef$factors$pvals) )
rownames(model_img_trait_dbrda_fit) <- c(names(model_img_trait_dbrda_ef$vectors$r),names(model_img_trait_dbrda_ef$factors$r))
model_img_trait_dbrda_fit
write.csv(model_img_trait_dbrda_fit, file="plots/img_trait_comp_dbrda_fit.csv", row.names=TRUE)

# Plot dbRDA
pdf(file="plots/img_trait_comp_dbrda_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_img_trait_dbrda_scores$sites[,1])-1, max(model_img_trait_dbrda_scores$sites[,1])+1),
	 ylim=c(min(model_img_trait_dbrda_scores$sites[,2]), max(model_img_trait_dbrda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_img_trait_dbrda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_img_trait_dbrda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_img_trait_dbrda, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_img_trait_dbrda, display="sites", labels=rownames(img_trait_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_img_trait_dbrda_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- Selected dbRDA of RGB x comp_list ----------
# Model with intercept only
model_0_img_trait_rda <- vegan::dbrda(formula=img_trait_list ~ 1, data=img_trait_rda_y, distance="bray")

# Model with all explanatory variables
model_1_img_trait_rda <- vegan::dbrda(formula=img_trait_list ~ ., data=img_trait_rda_y, distance="bray")

# Stepwise forward selection
#model_step_img_trait_rda <- ordistep(object=model_0_img_trait_rda, scope=formula(model_1_img_trait_rda), Pin=0.1, Pout=0.5, direction="forward", perm.max=100000)
model_step_img_trait_rda <- ordistep(object=model_0_img_trait_rda, scope=formula(model_1_img_trait_rda), Pin=0.05, Pout=0.1, trace=TRUE, direction="both", perm.max=100000, steps=10000)
model_step_img_trait_rda_scores <- vegan::scores(model_step_img_trait_rda)

# dbRDA with selected model by permutation tests in constrained ordination
model_img_trait_rda <- vegan::dbrda(formula=as.formula(model_step_img_trait_rda$terms), data=img_trait_rda_y, distance="bray")
model_img_trait_rda_ef_formula <- update(as.formula(model_step_img_trait_rda$terms), model_img_trait_rda ~ .)
model_img_trait_rda_ef_factors <- as.factor(sapply(strsplit(as.character(model_img_trait_rda_ef_formula)[[3]], "\\+"), function(x) { x <- gsub("(\\`|^ | $)","",x) }))
model_img_trait_rda_ef <- envfit(formula=model_img_trait_rda_ef_formula, data=img_trait_rda_y, perm=10000)
model_img_trait_rda_scores <- vegan::scores(model_img_trait_rda)

# Goodness of fit statistic: Squared correlation coefficient
model_img_trait_rda_fit <- data.frame(r2=c(model_img_trait_rda_ef$vectors$r,model_img_trait_rda_ef$factors$r),
									  pvals=c(model_img_trait_rda_ef$vectors$pvals,model_img_trait_rda_ef$factors$pvals) )
rownames(model_img_trait_rda_fit) <- c(names(model_img_trait_rda_ef$vectors$r),names(model_img_trait_rda_ef$factors$r))
model_img_trait_rda_fit
write.csv(model_img_trait_rda_fit, file="plots/img_trait_comp_dbrda_sel_fit.csv", row.names=TRUE)

# Plot results
pdf(file="plots/img_trait_comp_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_img_trait_rda_scores$sites[,1])-1, max(model_img_trait_rda_scores$sites[,1])+1),
	 ylim=c(min(model_img_trait_rda_scores$sites[,2]), max(model_img_trait_rda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_img_trait_rda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_img_trait_rda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_img_trait_rda, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_img_trait_rda, display="sites", labels=rownames(img_trait_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_img_trait_rda_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- Full dbRDA of RGB x class_list ----------
attach(as.data.frame(class_list))
img_trait_rda_formula <- formula(paste0("~ 0 + ", paste0('`',colnames(class_list),'`',collapse=" + ")))
img_trait_rda_y <- data.frame(model.matrix(img_trait_rda_formula))
detach(as.data.frame(class_list))

# Full dbRDA with compound list
model_img_trait_dbrda <- vegan::dbrda(formula=img_trait_list ~ ., data=img_trait_rda_y, distance="bray")
model_img_trait_dbrda_scores <- vegan::scores(model_img_trait_dbrda, choices=c(1,2))
model_img_trait_dbrda_ef <- envfit(model_img_trait_dbrda, class_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
model_img_trait_dbrda_fit <- data.frame(r2=c(model_img_trait_dbrda_ef$vectors$r,model_img_trait_dbrda_ef$factors$r),
										pvals=c(model_img_trait_dbrda_ef$vectors$pvals,model_img_trait_dbrda_ef$factors$pvals) )
rownames(model_img_trait_dbrda_fit) <- c(names(model_img_trait_dbrda_ef$vectors$r),names(model_img_trait_dbrda_ef$factors$r))
model_img_trait_dbrda_fit
write.csv(model_img_trait_dbrda_fit, file="plots/img_trait_class_dbrda_fit.csv", row.names=TRUE)

# Plot dbRDA
pdf(file="plots/img_trait_class_dbrda_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_img_trait_dbrda_scores$sites[,1])-1, max(model_img_trait_dbrda_scores$sites[,1])+1),
	 ylim=c(min(model_img_trait_dbrda_scores$sites[,2]), max(model_img_trait_dbrda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_img_trait_dbrda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_img_trait_dbrda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_img_trait_dbrda, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_img_trait_dbrda, display="sites", labels=rownames(img_trait_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_img_trait_dbrda_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- Selected dbRDA of RGB x class_list ----------
# Model with intercept only
model_0_img_trait_rda <- vegan::dbrda(formula=img_trait_list ~ 1, data=img_trait_rda_y, distance="bray")

# Model with all explanatory variables
model_1_img_trait_rda <- vegan::dbrda(formula=img_trait_list ~ ., data=img_trait_rda_y, distance="bray")

# Stepwise forward selection
model_step_img_trait_rda <- ordistep(object=model_0_img_trait_rda, scope=formula(model_1_img_trait_rda), Pin=0.05, Pout=0.1, trace=TRUE, direction="both", perm.max=100000, steps=10000)
model_step_img_trait_rda_scores <- vegan::scores(model_step_img_trait_rda)

# dbRDA with selected model by permutation tests in constrained ordination
model_img_trait_rda <- vegan::dbrda(formula=as.formula(model_step_img_trait_rda$terms), data=img_trait_rda_y, distance="bray")
model_img_trait_rda_ef_formula <- update(as.formula(model_step_img_trait_rda$terms), model_img_trait_rda ~ .)
model_img_trait_rda_ef_factors <- as.factor(sapply(strsplit(as.character(model_img_trait_rda_ef_formula)[[3]], "\\+"), function(x) { x <- gsub("(\\`|^ | $)","",x) }))
model_img_trait_rda_ef <- envfit(formula=model_img_trait_rda_ef_formula, data=img_trait_rda_y, perm=10000)
model_img_trait_rda_scores <- vegan::scores(model_img_trait_rda)

# Goodness of fit statistic: Squared correlation coefficient
model_img_trait_rda_fit <- data.frame(r2=c(model_img_trait_rda_ef$vectors$r,model_img_trait_rda_ef$factors$r),
									  pvals=c(model_img_trait_rda_ef$vectors$pvals,model_img_trait_rda_ef$factors$pvals) )
rownames(model_img_trait_rda_fit) <- c(names(model_img_trait_rda_ef$vectors$r),names(model_img_trait_rda_ef$factors$r))
model_img_trait_rda_fit
write.csv(model_img_trait_rda_fit, file="plots/img_trait_class_dbrda_sel_fit.csv", row.names=TRUE)

# Plot results
pdf(file="plots/img_trait_class_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_img_trait_rda_scores$sites[,1])-1, max(model_img_trait_rda_scores$sites[,1])+1),
	 ylim=c(min(model_img_trait_rda_scores$sites[,2]), max(model_img_trait_rda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_img_trait_rda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_img_trait_rda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_img_trait_rda, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_img_trait_rda, display="sites", labels=rownames(img_trait_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_img_trait_rda_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- Full dbRDA of RGB x superclass_list ----------
attach(as.data.frame(superclass_list))
img_trait_rda_formula <- formula(paste0("~ 0 + ", paste0('`',colnames(superclass_list),'`',collapse=" + ")))
img_trait_rda_y <- data.frame(model.matrix(img_trait_rda_formula))
detach(as.data.frame(superclass_list))

# Full dbRDA with compound list
model_img_trait_dbrda <- vegan::dbrda(formula=img_trait_list ~ ., data=img_trait_rda_y, distance="bray")
model_img_trait_dbrda_scores <- vegan::scores(model_img_trait_dbrda, choices=c(1,2))
model_img_trait_dbrda_ef <- envfit(model_img_trait_dbrda, superclass_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
model_img_trait_dbrda_fit <- data.frame(r2=c(model_img_trait_dbrda_ef$vectors$r,model_img_trait_dbrda_ef$factors$r),
										pvals=c(model_img_trait_dbrda_ef$vectors$pvals,model_img_trait_dbrda_ef$factors$pvals) )
rownames(model_img_trait_dbrda_fit) <- c(names(model_img_trait_dbrda_ef$vectors$r),names(model_img_trait_dbrda_ef$factors$r))
model_img_trait_dbrda_fit
write.csv(model_img_trait_dbrda_fit, file="plots/img_trait_superclass_dbrda_fit.csv", row.names=TRUE)

# Plot dbRDA
pdf(file="plots/img_trait_superclass_dbrda_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_img_trait_dbrda_scores$sites[,1])-1, max(model_img_trait_dbrda_scores$sites[,1])+1),
	 ylim=c(min(model_img_trait_dbrda_scores$sites[,2]), max(model_img_trait_dbrda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_img_trait_dbrda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_img_trait_dbrda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_img_trait_dbrda, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_img_trait_dbrda, display="sites", labels=rownames(img_trait_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_img_trait_dbrda_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- Selected dbRDA of RGB x superclass_list ----------
# Model with intercept only
model_0_img_trait_rda <- vegan::dbrda(formula=img_trait_list ~ 1, data=img_trait_rda_y, distance="bray")

# Model with all explanatory variables
model_1_img_trait_rda <- vegan::dbrda(formula=img_trait_list ~ ., data=img_trait_rda_y, distance="bray")

# Stepwise forward selection
model_step_img_trait_rda <- ordistep(object=model_0_img_trait_rda, scope=formula(model_1_img_trait_rda), Pin=0.5, Pout=0.7, trace=TRUE, direction="both", perm.max=100000, steps=10000)
model_step_img_trait_rda_scores <- vegan::scores(model_step_img_trait_rda)

# dbRDA with selected model by permutation tests in constrained ordination
model_img_trait_rda <- vegan::dbrda(formula=as.formula(model_step_img_trait_rda$terms), data=img_trait_rda_y, distance="bray")
model_img_trait_rda_ef_formula <- update(as.formula(model_step_img_trait_rda$terms), model_img_trait_rda ~ .)
model_img_trait_rda_ef_factors <- as.factor(sapply(strsplit(as.character(model_img_trait_rda_ef_formula)[[3]], "\\+"), function(x) { x <- gsub("(\\`|^ | $)","",x) }))
model_img_trait_rda_ef <- envfit(formula=model_img_trait_rda_ef_formula, data=img_trait_rda_y, perm=10000)
model_img_trait_rda_scores <- vegan::scores(model_img_trait_rda)

# Goodness of fit statistic: Squared correlation coefficient
model_img_trait_rda_fit <- data.frame(r2=c(model_img_trait_rda_ef$vectors$r,model_img_trait_rda_ef$factors$r),
									  pvals=c(model_img_trait_rda_ef$vectors$pvals,model_img_trait_rda_ef$factors$pvals) )
rownames(model_img_trait_rda_fit) <- c(names(model_img_trait_rda_ef$vectors$r),names(model_img_trait_rda_ef$factors$r))
model_img_trait_rda_fit
write.csv(model_img_trait_rda_fit, file="plots/img_trait_superclass_dbrda_sel_fit.csv", row.names=TRUE)

# Plot results
pdf(file="plots/img_trait_superclass_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(model_img_trait_rda_scores$sites[,1])-1, max(model_img_trait_rda_scores$sites[,1])+1),
	 ylim=c(min(model_img_trait_rda_scores$sites[,2]), max(model_img_trait_rda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(model_img_trait_rda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(model_img_trait_rda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(model_img_trait_rda, display="sites", pch=19, col=mzml_pheno_colors_samples)
text(model_img_trait_rda, display="sites", labels=rownames(img_trait_list), col=mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(model_img_trait_rda_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ############################## Chemotaxonomy Ingroup ##############################



# ---------- Taxonomic tree from genetic data ----------
# Read phylogenetic tree
phylo_ingroup_tree <- read.tree("data/trnLF_ingroup_newick.txt")
phylo_ingroup_tree <- chronos(phylo_ingroup_tree)
phylo_ingroup_heights <- phylo_ingroup_tree$edge.length
phylo_ingroup_dend <- as.dendrogram(phylo_ingroup_tree)
labels(phylo_ingroup_dend) <- c("R.sorocarpa", "R.glauca", "R.warnstorfii")

#plot(as.phylo(phylo_ingroup_dend), type="phylogram", direction="rightwards", align.tip.label=T, use.edge.length=T, show.tip.label=T, show.node.label=TRUE, add.scale.bar=TRUE)
#plot(as.hclust(phylo_ingroup_dend), show.node.label=TRUE, add.scale.bar=TRUE, hang=-1)
#with(pvclust:::hc2axes(as.hclust(phylo_ingroup_dend)), text(x.axis, y.axis, y.axis, pos=1))
#ggtree(as.phylo(phylo_ingroup_dend)) + geom_tiplab(size=4) + geom_text2(aes(subset=!isTip, label=node), size=5, color="black", hjust=1, vjust=1) 

# Distance matrix of phylogenetic tree using Bray-Curtis
phylo_ingroup_dist <- cophenetic.phylo(phylo_ingroup_tree)
phylo_ingroup_dist <- vegdist(phylo_ingroup_dist, method="bray")

# Hierarchical clustering
phylo_ingroup_hclust <- hclust(phylo_ingroup_dist, method="complete")



# ---------- Chemotaxonomic tree for feat_list ----------
# Merge feat_list for species from samples
phylo_feat_list <- NULL
for (i in unique(mzml_pheno_samples)) phylo_feat_list <- rbind(phylo_feat_list, apply(X=feat_list[mzml_pheno_samples==i,], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_feat_list) <- unique(mzml_pheno_samples)

# Distance matrix of feat_list using Bray-Curtis
phylo_feat_dist <- vegdist(phylo_feat_list, method="bray")

# Hierarchical clustering
phylo_feat_hclust <- hclust(phylo_feat_dist, method="complete")

# Optimal order
phylo_feat_opti <- order.optimal(phylo_feat_dist, phylo_feat_hclust$merge)
phylo_feat_oclust <- phylo_feat_hclust
phylo_feat_oclust$merge <- phylo_feat_opti$merge
phylo_feat_oclust$order <- phylo_feat_opti$order

# Dendrogram
phylo_feat_dend <- as.dendrogram(phylo_feat_oclust)



# ---------- Chemotaxonomic tree for compound list ----------
# Merge comp_list for species from samples
phylo_comp_list <- NULL
for (i in unique(mzml_pheno_samples)) phylo_comp_list <- rbind(phylo_comp_list, apply(X=comp_list[mzml_pheno_samples==i,], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_comp_list) <- unique(mzml_pheno_samples)

# Distance matrix of feat_list using Bray-Curtis
phylo_comp_dist <- vegdist(phylo_comp_list, method="bray")

# Hierarchical clustering
phylo_comp_hclust <- hclust(phylo_comp_dist, method="complete")

# Optimal order
phylo_comp_opti <- order.optimal(phylo_comp_dist, phylo_comp_hclust$merge)
phylo_comp_oclust <- phylo_comp_hclust
phylo_comp_oclust$merge <- phylo_comp_opti$merge
phylo_comp_oclust$order <- phylo_comp_opti$order

# Dendrogram
phylo_comp_dend <- as.dendrogram(phylo_comp_oclust)



# ---------- Chemotaxonomic tree for class list ----------
# Merge class_list for species from samples
phylo_class_list <- NULL
for (i in unique(mzml_pheno_samples)) phylo_class_list <- rbind(phylo_class_list, apply(X=class_list[mzml_pheno_samples==i,], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_class_list) <- unique(mzml_pheno_samples)

# Distance matrix of feat_list using Bray-Curtis
phylo_class_dist <- vegdist(phylo_class_list, method="bray")

# Hierarchical clustering
phylo_class_hclust <- hclust(phylo_class_dist, method="complete")

# Optimal order
phylo_class_opti <- order.optimal(phylo_class_dist, phylo_class_hclust$merge)
phylo_class_oclust <- phylo_class_hclust
phylo_class_oclust$merge <- phylo_class_opti$merge
phylo_class_oclust$order <- phylo_class_opti$order

# Dendrogram
phylo_class_dend <- as.dendrogram(phylo_class_oclust)

# Reorder dendrogram
#temp <- phylo_class_dend[[1]]
#phylo_class_dend[[1]] <- phylo_class_dend[[2]]
#phylo_class_dend[[2]] <- temp
#rm(temp)



# ---------- Chemotaxonomic tree for superclass list ----------
# Merge superclass_list for species from samples
phylo_superclass_list <- NULL
for (i in unique(mzml_pheno_samples)) phylo_superclass_list <- rbind(phylo_superclass_list, apply(X=superclass_list[mzml_pheno_samples==i,], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_superclass_list) <- unique(mzml_pheno_samples)

# Distance matrix of feat_list using Bray-Curtis
phylo_superclass_dist <- vegdist(phylo_superclass_list, method="bray")

# Hierarchical clustering
phylo_superclass_hclust <- hclust(phylo_superclass_dist, method="complete")

# Optimal order
phylo_superclass_opti <- order.optimal(phylo_superclass_dist, phylo_superclass_hclust$merge)
phylo_superclass_oclust <- phylo_superclass_hclust
phylo_superclass_oclust$merge <- phylo_superclass_opti$merge
phylo_superclass_oclust$order <- phylo_superclass_opti$order

# Dendrogram
phylo_superclass_dend <- as.dendrogram(phylo_superclass_oclust)



# ---------- Chemotaxonomic tree for NP class list ----------
# Merge npclass_list for species from samples
phylo_npclass_list <- NULL
for (i in unique(mzml_pheno_samples)) phylo_npclass_list <- rbind(phylo_npclass_list, apply(X=npclass_list[mzml_pheno_samples==i,], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_npclass_list) <- unique(mzml_pheno_samples)

# Distance matrix of feat_list using Bray-Curtis
phylo_npclass_dist <- vegdist(phylo_npclass_list, method="bray")

# Hierarchical clustering
phylo_npclass_hclust <- hclust(phylo_npclass_dist, method="complete")

# Optimal order
phylo_npclass_opti <- order.optimal(phylo_npclass_dist, phylo_npclass_hclust$merge)
phylo_npclass_oclust <- phylo_npclass_hclust
phylo_npclass_oclust$merge <- phylo_npclass_opti$merge
phylo_npclass_oclust$order <- phylo_npclass_opti$order

# Dendrogram
phylo_npclass_dend <- as.dendrogram(phylo_npclass_oclust)



# ---------- Chemotaxonomic tree for molecular descriptors list ----------
# Merge mdes_list for species from samples
phylo_mdes_list <- NULL
for (i in unique(mzml_pheno_samples)) phylo_mdes_list <- rbind(phylo_mdes_list, apply(X=mdes_list[mzml_pheno_samples==i,], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_mdes_list) <- unique(mzml_pheno_samples)

# Distance matrix of feat_list using Bray-Curtis
phylo_mdes_dist <- vegdist(phylo_mdes_list, method="bray")

# Hierarchical clustering
phylo_mdes_hclust <- hclust(phylo_mdes_dist, method="complete")

# Optimal order
phylo_mdes_opti <- order.optimal(phylo_mdes_dist, phylo_mdes_hclust$merge)
phylo_mdes_oclust <- phylo_mdes_hclust
phylo_mdes_oclust$merge <- phylo_mdes_opti$merge
phylo_mdes_oclust$order <- phylo_mdes_opti$order

# Dendrogram
phylo_mdes_dend <- as.dendrogram(phylo_mdes_oclust)



# ---------- Plot chemotaxonomic trees ----------
# Plot chemotaxonomic trees
pdf("plots/chemotax_phylo_ingroup_trees.pdf", encoding="ISOLatin1", pointsize=8, width=5*2, height=2.8, family="Helvetica")
par(mfrow=c(1,5), mar=c(1,1,2,1), cex=1.0)
phylo_max <- 3.00; plot(as.phylo(phylo_ingroup_dend), type="phylogram", direction="leftwards", x.lim=c(0,phylo_max), label.offset=phylo_max/30, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=c(mzml_pheno_colors[2],mzml_pheno_colors[1],mzml_pheno_colors[3]), font=4, main="")
edgelabels(text=round(node.depth.edgelength(as.phylo(phylo_ingroup_dend))[-4],3), col="black", bg="white", frame="none", adj=c(0.5,-0.5), cex=0.8)
mtext(text="(a)", adj=0, line=0.5, font=2, cex=1.2)
phylo_max <- max(phylo_comp_dist); plot(as.phylo(phylo_comp_oclust), type="phylogram", direction="rightwards", x.lim=c(0,phylo_max), label.offset=phylo_max/30, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=mzml_pheno_colors, font=4, main="")
edgelabels(text=round(node.depth.edgelength(as.phylo(phylo_comp_oclust))[-4],3), col="black", bg="white", frame="none", adj=c(0.5,-0.5), cex=0.8)
mtext(text="(b)", adj=0, line=0.5, font=2, cex=1.2)
phylo_max <- max(phylo_class_dist); plot(as.phylo(phylo_class_oclust), type="phylogram", direction="rightwards", x.lim=c(0,phylo_max), label.offset=phylo_max/30, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=mzml_pheno_colors, font=4, main="")
edgelabels(text=round(node.depth.edgelength(as.phylo(phylo_class_oclust))[-4],3), col="black", bg="white", frame="none", adj=c(0.5,-0.5), cex=0.8)
mtext(text="(c)", adj=0, line=0.5, font=2, cex=1.2)
phylo_max <- max(phylo_superclass_dist); plot(as.phylo(phylo_superclass_oclust), type="phylogram", direction="rightwards", x.lim=c(0,phylo_max), label.offset=phylo_max/30, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=mzml_pheno_colors, font=4, main="")
edgelabels(text=round(node.depth.edgelength(as.phylo(phylo_superclass_oclust))[-4],3), col="black", bg="white", frame="none", adj=c(0.5,-0.5), cex=0.8)
mtext(text="(d)", adj=0, line=0.5, font=2, cex=1.2)
phylo_max <- max(phylo_mdes_dist); plot(as.phylo(phylo_mdes_oclust), type="phylogram", direction="rightwards", x.lim=c(0,phylo_max), label.offset=phylo_max/30, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=mzml_pheno_colors, font=4, main="")
edgelabels(text=round(node.depth.edgelength(as.phylo(phylo_mdes_oclust))[-4],3), col="black", bg="white", frame="none", adj=c(0.5,-0.5), cex=0.8)
mtext(text="(e)", adj=0, line=0.5, font=2, cex=1.2)
dev.off()



# ---------- Clustering metrics ----------
phylo_model_tree_metrics <- NULL
phylo_model_tree_metrics <- rbind(phylo_model_tree_metrics, c("mantel"=round(mantel(xdis=phylo_ingroup_dist, ydis=phylo_comp_dist, method="pearson", permutations=10000)$statistic, 3),
															  "cor"=cor(phylo_ingroup_dist, phylo_comp_dist, method="pearson"),
															  "cop"=cor_cophenetic(hclust(phylo_ingroup_dist), hclust(phylo_comp_dist), method="pearson"),
															  "robinson-foulds"=RF.dist(as.phylo(phylo_ingroup_dend), as.phylo(phylo_comp_oclust), normalize=TRUE, check.labels=FALSE, rooted=TRUE)
))
phylo_model_tree_metrics <- rbind(phylo_model_tree_metrics, c("mantel"=round(mantel(xdis=phylo_ingroup_dist, ydis=phylo_class_dist, method="pearson", permutations=10000)$statistic, 3),
															  "cor"=cor(phylo_ingroup_dist, phylo_class_dist, method="pearson"),
															  "cop"=cor_cophenetic(hclust(phylo_ingroup_dist), hclust(phylo_class_dist), method="pearson"),
															  "robinson-foulds"=RF.dist(as.phylo(phylo_ingroup_dend), as.phylo(phylo_class_oclust), normalize=TRUE, check.labels=FALSE, rooted=TRUE)
))
phylo_model_tree_metrics <- rbind(phylo_model_tree_metrics, c("mantel"=round(mantel(xdis=phylo_ingroup_dist, ydis=phylo_superclass_dist, method="pearson", permutations=10000)$statistic, 3),
															  "cor"=cor(phylo_ingroup_dist, phylo_superclass_dist, method="pearson"),
															  "cop"=cor_cophenetic(hclust(phylo_ingroup_dist), hclust(phylo_superclass_dist), method="pearson"),
															  "robinson-foulds"=RF.dist(as.phylo(phylo_ingroup_dend), as.phylo(phylo_superclass_oclust), normalize=TRUE, check.labels=FALSE, rooted=TRUE)
))
phylo_model_tree_metrics <- rbind(phylo_model_tree_metrics, c("mantel"=round(mantel(xdis=phylo_ingroup_dist, ydis=phylo_mdes_dist, method="pearson", permutations=10000)$statistic, 3),
															  "cor"=cor(phylo_ingroup_dist, phylo_mdes_dist, method="pearson"),
															  "cop"=cor_cophenetic(hclust(phylo_ingroup_dist), hclust(phylo_mdes_dist), method="pearson"),
															  "robinson-foulds"=RF.dist(as.phylo(phylo_ingroup_dend), as.phylo(phylo_mdes_oclust), normalize=TRUE, check.labels=FALSE, rooted=TRUE)
))
rownames(phylo_model_tree_metrics) <- c("comp", "class", "superclass", "mdes")
write.csv(phylo_model_tree_metrics, file="plots/chemotax_phylo_ingroup_trees.csv", row.names=TRUE)



# ############################## Chemotaxonomy with outgroup ##############################



# ---------- Peak detection ----------
source("phylo_peak_detection_neg.r")
source("phylo_peak_detection_pos.r")



# ---------- Create merged pos+neg objects ----------
# factors
phylo_mzml_names <- gsub(x=phylo_mzml_names_pos, pattern="\\.pos.*", replacement="")
phylo_mzml_pheno_samples <- gsub(x=phylo_mzml_pheno_samples_pos, pattern="R\\..*", replacement="Riccia.spp.")
phylo_mzml_pheno_colors <- c("forestgreen", "seagreen3")
phylo_mzml_pheno_colors_samples <- c("forestgreen", "forestgreen", "forestgreen", "seagreen3", "seagreen3", "seagreen3", "seagreen3", "seagreen3", "seagreen3", "seagreen3", "seagreen3", "seagreen3")

# phylo_feat_list
phylo_feat_list <- cbind(phylo_feat_list_pos, phylo_feat_list_neg)
rownames(phylo_feat_list) <- gsub(x=rownames(phylo_feat_list_pos), pattern="\\.pos.*", replacement="")
colnames(phylo_feat_list) <- c(paste0(colnames(phylo_feat_list_pos),"_pos"), paste0(colnames(phylo_feat_list_neg),"_neg"))

# phylo_comp_list
phylo_comp_list <- cbind(phylo_feat_list_pos[, c(rownames(phylo_ms1_def_pos)[phylo_ms1_def_pos$has_ms2==1])],
						 phylo_feat_list_neg[, c(rownames(phylo_ms1_def_neg)[phylo_ms1_def_neg$has_ms2==1])] )
colnames(phylo_comp_list) <- c(paste0(rownames(phylo_ms1_def_pos)[phylo_ms1_def_pos$has_ms2==1], "_pos"),
							   paste0(rownames(phylo_ms1_def_neg)[phylo_ms1_def_neg$has_ms2==1], "_neg") )

# phylo_bina_list
phylo_bina_list <- cbind(phylo_bina_list_pos[, c(rownames(phylo_ms1_def_pos)[phylo_ms1_def_pos$has_ms2==1])],
						 phylo_bina_list_neg[, c(rownames(phylo_ms1_def_neg)[phylo_ms1_def_neg$has_ms2==1])] )
rownames(phylo_bina_list) <- gsub(x=rownames(phylo_bina_list_pos), pattern="\\.pos.*", replacement="")
colnames(phylo_bina_list) <- c(paste0(rownames(phylo_ms1_def_pos)[phylo_ms1_def_pos$has_ms2==1], "_pos"),
							   paste0(rownames(phylo_ms1_def_neg)[phylo_ms1_def_neg$has_ms2==1], "_neg") )

# phylo_uniq_list
phylo_uniq_list <- cbind(phylo_uniq_list_pos[, c(rownames(phylo_ms1_def_pos)[phylo_ms1_def_pos$has_ms2==1])],
						 phylo_uniq_list_neg[, c(rownames(phylo_ms1_def_neg)[phylo_ms1_def_neg$has_ms2==1])] )
rownames(phylo_uniq_list) <- gsub(x=rownames(phylo_uniq_list_pos), pattern="\\.pos.*", replacement="")
colnames(phylo_uniq_list) <- c(paste0(rownames(phylo_ms1_def_pos)[phylo_ms1_def_pos$has_ms2==1], "_pos"),
							   paste0(rownames(phylo_ms1_def_neg)[phylo_ms1_def_neg$has_ms2==1], "_neg") )

# Create data frame
phylo_model_div             <- data.frame(features=apply(X=phylo_bina_list, MARGIN=1, FUN=function(x) { sum(x) } ))
phylo_model_div$richness    <- apply(X=phylo_bina_list, MARGIN=1, FUN=function(x) { sum(x) } )
#phylo_model_div$menhinick   <- apply(X=phylo_bina_list, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
phylo_model_div$shannon     <- apply(X=phylo_comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
phylo_model_div$pielou      <- apply(X=scale(phylo_comp_list, center=F, scale=T), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#phylo_model_div$chao        <- vegan::specpool2vect(X=vegan::specpool(phylo_feat_list, species), index="chao")
phylo_model_div$simpson     <- apply(X=phylo_comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
phylo_model_div$inverse     <- apply(X=phylo_comp_list, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
phylo_model_div$fisher      <- apply(X=phylo_comp_list, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
phylo_model_div$unique      <- apply(X=phylo_uniq_list, MARGIN=1, FUN=function(x) { sum(x) })

# Remove NAs if present
phylo_model_div[is.na(phylo_model_div)] <- 0

# phylo_class_list
phylo_class_list <- cbind(phylo_class_list_pos, phylo_class_list_neg)
rownames(phylo_class_list) <- gsub(x=rownames(phylo_class_list_pos), pattern="\\.pos.*", replacement="")
phylo_class_list <- cbind(sapply(unique(colnames(phylo_class_list)[duplicated(colnames(phylo_class_list))]), function(x) rowSums(phylo_class_list[,grepl(paste(x, "$", sep=""), colnames(phylo_class_list))])), phylo_class_list[,!duplicated(colnames(phylo_class_list)) & ! duplicated(colnames(phylo_class_list), fromLast=TRUE)])

# phylo_class_int_list
phylo_class_int_list <- cbind(phylo_class_int_list_pos, phylo_class_int_list_neg)
rownames(phylo_class_int_list) <- gsub(x=rownames(phylo_class_int_list_pos), pattern="\\.pos.*", replacement="")
phylo_class_int_list <- cbind(sapply(unique(colnames(phylo_class_int_list)[duplicated(colnames(phylo_class_int_list))]), function(x) rowSums(phylo_class_int_list[,grepl(paste(x, "$", sep=""), colnames(phylo_class_int_list))])), phylo_class_int_list[,!duplicated(colnames(phylo_class_int_list)) & ! duplicated(colnames(phylo_class_int_list), fromLast=TRUE)])

# phylo_superclass_list
phylo_superclass_list <- cbind(phylo_superclass_list_pos, phylo_superclass_list_neg)
rownames(phylo_superclass_list) <- gsub(x=rownames(phylo_superclass_list_pos), pattern="\\.pos.*", replacement="")
phylo_superclass_list <- cbind(sapply(unique(colnames(phylo_superclass_list)[duplicated(colnames(phylo_superclass_list))]), function(x) rowSums(phylo_superclass_list[,grepl(paste(x, "$", sep=""), colnames(phylo_superclass_list))])), phylo_superclass_list[,!duplicated(colnames(phylo_superclass_list)) & ! duplicated(colnames(phylo_superclass_list), fromLast=TRUE)])

# phylo_superclass_int_list
phylo_superclass_int_list <- cbind(phylo_superclass_int_list_pos, phylo_superclass_int_list_neg)
rownames(phylo_superclass_int_list) <- gsub(x=rownames(phylo_superclass_int_list_pos), pattern="\\.pos.*", replacement="")
phylo_superclass_int_list <- cbind(sapply(unique(colnames(phylo_superclass_int_list)[duplicated(colnames(phylo_superclass_int_list))]), function(x) rowSums(phylo_superclass_int_list[,grepl(paste(x, "$", sep=""), colnames(phylo_superclass_int_list))])), phylo_superclass_int_list[,!duplicated(colnames(phylo_superclass_int_list)) & ! duplicated(colnames(phylo_superclass_int_list), fromLast=TRUE)])



# ---------- Histogram ----------
pdf(file="plots/phylo_ms1_phylo_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(phylo_feat_list))
dev.off()



# ---------- Variation partitioning ----------
model_varpart <- varpart(scale(phylo_feat_list), ~ phylo_mzml_pheno_samples, ~ phylo_mzml_pheno_samples)

# Plot results
pdf(file="plots/phylo_ms1_phylo_feat_list_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca <- prcomp(phylo_feat_list, scale=TRUE, center=TRUE)

# PC-Axes 1+2
pdf(paste("plots/phylo_ms1_phylo_feat_list_pca12.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 1], model_pca$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC1: ", format(summary(model_pca)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples, bg=phylo_mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 1], model_pca$x[,2], labels=phylo_mzml_names, col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()

# PC-Axes 2+3
pdf(paste("plots/phylo_ms1_phylo_feat_list_pca23.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 2], model_pca$x[,3], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples, bg=phylo_mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 2], model_pca$x[,3], labels=phylo_mzml_names, col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()



# ---------- Diversity ----------
# Plot unique features
pdf(paste("plots/phylo_ms1_phylo_feat_list_richness_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$unique ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Number of unique compounds", xlab="treatment", ylab="number of unique compounds")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$unique, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot total features
pdf(paste("plots/phylo_ms1_phylo_feat_list_richness_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$features ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Number of compounds", xlab="treatment", ylab="number of compounds")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$features, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Shannon index
pdf(paste("plots/phylo_ms1_phylo_feat_list_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$shannon ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$shannon, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Pielou evenness for species
pdf(paste("plots/phylo_ms1_phylo_feat_list_diversity_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$pielou ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Pielou\'s evenness", xlab="treatment", ylab="Pielou diversity index (J)")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$pielou, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()



# ---------- Variable Selection ----------
# PLS
phylo_sel_pls <- f.select_features_pls(feat_matrix=phylo_feat_list, sel_factor=phylo_mzml_pheno_samples, sel_colors=phylo_mzml_pheno_colors, components=2, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/phylo_ms1_phylo_feat_list_select_pls_roc.pdf")
print(paste("Number of selected compounds:", f.count.selected_features(sel_feat=phylo_sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=phylo_feat_list, sel_feat=phylo_sel_pls$`_selected_variables_`, sample_colors=phylo_mzml_pheno_colors_samples, filename="plots/phylo_ms1_phylo_feat_list_select_pls.pdf", main="PLS")
phylo_sel_pls$`_selected_variables_`
phylo_sel_pls$`_model_r2_`
phylo_sel_pls$`_multiclass_metrics_`



# ---------- Histogram ----------
pdf(file="plots/phylo_ms1_phylo_comp_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(phylo_comp_list))
dev.off()



# ---------- Variation partitioning ----------
model_varpart <- varpart(scale(phylo_comp_list), ~ phylo_mzml_pheno_samples, ~ phylo_mzml_pheno_samples)

# Plot results
pdf(file="plots/phylo_ms1_phylo_comp_list_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca <- prcomp(phylo_comp_list, scale=TRUE, center=TRUE)

# PC-Axes 1+2
pdf(paste("plots/phylo_ms1_phylo_comp_list_pca12.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 1], model_pca$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC1: ", format(summary(model_pca)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples, bg=phylo_mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 1], model_pca$x[,2], labels=phylo_mzml_names, col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()

# PC-Axes 2+3
pdf(paste("plots/phylo_ms1_phylo_comp_list_pca23.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 2], model_pca$x[,3], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples, bg=phylo_mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 2], model_pca$x[,3], labels=phylo_mzml_names, col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()



# ---------- Diversity ----------
# Plot unique features
pdf(paste("plots/phylo_ms1_phylo_comp_list_richness_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$unique ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Number of unique compounds", xlab="treatment", ylab="number of unique compounds")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$unique, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot total features
pdf(paste("plots/phylo_ms1_phylo_comp_list_richness_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$features ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Number of compounds", xlab="treatment", ylab="number of compounds")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$features, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Shannon index
pdf(paste("plots/phylo_ms1_phylo_comp_list_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$shannon ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$shannon, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Pielou evenness for species
pdf(paste("plots/phylo_ms1_phylo_comp_list_diversity_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$pielou ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Pielou\'s evenness", xlab="treatment", ylab="Pielou diversity index (J)")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$pielou, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()



# ---------- Variable Selection ----------
# Load identification
riccia_id_outgroup <- read.table(file="data/riccia_id_outgroup.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)
riccia_id_outgroup[is.na(riccia_id_outgroup)] <- "-"
riccia_id_outgroup$sel_names <- riccia_id_outgroup$XCMS.name 

for (i in 1:length(riccia_id_outgroup$sel_names)) {
	if ((riccia_id_outgroup$Compound.Name[i] != "-") & (riccia_id_outgroup$Compound.Name[i] != "n.d") & (riccia_id_outgroup$Compound.Name[i] != "")) {
		riccia_id_outgroup$sel_names[i] <- substr(riccia_id_outgroup$Compound.Name[i], 1, 50)
	}
}

# PLS
phylo_sel_pls <- f.select_features_pls(feat_matrix=phylo_comp_list, sel_factor=phylo_mzml_pheno_samples, sel_colors=phylo_mzml_pheno_colors, components=2, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/phylo_ms1_phylo_comp_list_select_pls_roc.pdf")
print(paste("Number of selected compounds:", f.count.selected_features(sel_feat=phylo_sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=phylo_comp_list, sel_feat=phylo_sel_pls$`_selected_variables_`, sel_names=riccia_id_outgroup$sel_names, sample_colors=phylo_mzml_pheno_colors_samples, plot_width=14, plot_height=5, filename="plots/phylo_ms1_phylo_comp_list_select_pls.pdf", main="")
f.heatmap.selected_features(feat_list=phylo_comp_list, sel_feat=phylo_sel_pls$`_selected_variables_`, sel_names=paste0("     ",phylo_sel_pls$`_selected_variables_`), sample_colors=phylo_mzml_pheno_colors_samples, plot_width=14, plot_height=4, filename="plots/phylo_ms1_phylo_comp_list_select_pls.pdf", main="")
phylo_sel_pls$`_selected_variables_`
phylo_sel_pls$`_model_r2_`
phylo_sel_pls$`_multiclass_metrics_`
#boxplot(phylo_comp_list[,"FT0024_pos"] ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors)



# ---------- Histogram ----------
pdf(file="plots/phylo_ms1_phylo_class_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(as.matrix(phylo_class_list)))
dev.off()



# ---------- Variation partitioning ----------
model_varpart <- varpart(scale(phylo_class_list), ~ phylo_mzml_pheno_samples, ~ phylo_mzml_pheno_samples)

# Plot results
pdf(file="plots/phylo_ms1_phylo_class_list_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca <- prcomp(phylo_class_list, scale=TRUE, center=TRUE)

# PC-Axes 1+2
pdf(paste("plots/phylo_ms1_phylo_class_list_pca12.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 1], model_pca$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC1: ", format(summary(model_pca)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples, bg=phylo_mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 1], model_pca$x[,2], labels=phylo_mzml_names, col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()

# PC-Axes 2+3
pdf(paste("plots/phylo_ms1_phylo_class_list_pca23.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 2], model_pca$x[,3], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples, bg=phylo_mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 2], model_pca$x[,3], labels=phylo_mzml_names, col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()



# ---------- Diversity ----------
# Plot unique features
pdf(paste("plots/phylo_ms1_phylo_class_list_richness_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$unique ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Number of unique compounds", xlab="treatment", ylab="number of unique compounds")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$unique, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot total features
pdf(paste("plots/phylo_ms1_phylo_class_list_richness_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$features ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Number of compounds", xlab="treatment", ylab="number of compounds")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$features, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Shannon index
pdf(paste("plots/phylo_ms1_phylo_class_list_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$shannon ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$shannon, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Pielou evenness for species
pdf(paste("plots/phylo_ms1_phylo_class_list_diversity_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(phylo_model_div$pielou ~ phylo_mzml_pheno_samples, col=phylo_mzml_pheno_colors, names=NA, main="Pielou\'s evenness", xlab="treatment", ylab="Pielou diversity index (J)")
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(phylo_mzml_pheno_samples), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=phylo_model_div$pielou, term=as.factor(phylo_mzml_pheno_samples))
text(1:length(unique(phylo_mzml_pheno_samples)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()



# ---------- Variable Selection ----------
# PLS
phylo_sel_pls <- f.select_features_pls(feat_matrix=phylo_class_list, sel_factor=phylo_mzml_pheno_samples, sel_colors=phylo_mzml_pheno_colors, components=2, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/phylo_ms1_phylo_class_list_select_pls_roc.pdf")
print(paste("Number of selected compounds:", f.count.selected_features(sel_feat=phylo_sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=phylo_class_list, sel_feat=phylo_sel_pls$`_selected_variables_`, sample_colors=phylo_mzml_pheno_colors_samples, plot_width=3, plot_height=5, filename="plots/phylo_ms1_phylo_class_list_select_pls.pdf", main="")
phylo_sel_pls$`_selected_variables_`
phylo_sel_pls$`_model_r2_`
phylo_sel_pls$`_multiclass_metrics_`



# ---------- Sunburst plot of classes of all samples ----------
# Riccia
phylo_riccia_div_classes <- merge(phylo_riccia_div_classes_pos, phylo_riccia_div_classes_neg, by="row.names", all.x=TRUE, all.y=TRUE)
colnames(phylo_riccia_div_classes) <- c("superclass", "phylo_riccia_div_classes_pos", "phylo_riccia_div_classes_neg")
phylo_riccia_div_classes[is.na(phylo_riccia_div_classes)] <- 0
phylo_riccia_div_classes$frequency <- apply(X=phylo_riccia_div_classes[,c(2:ncol(phylo_riccia_div_classes))], MARGIN=1, FUN=function(x) { sum(x) })
rownames(phylo_riccia_div_classes) <- phylo_riccia_div_classes[,1]
phylo_riccia_div_classes <- phylo_riccia_div_classes[-which(phylo_riccia_div_classes$frequency==0),]

pdf(file="plots/phylo_riccia_ms2_classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(phylo_riccia_div_classes), phylo_riccia_div_classes$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(phylo_riccia_div_classes), "Number of spectra"=phylo_riccia_div_classes$frequency), file="plots/phylo_riccia_ms2_classes_sunburst.csv", row.names=FALSE)

# Lunularia
phylo_lunularia_div_classes <- merge(phylo_lunularia_div_classes_pos, phylo_lunularia_div_classes_neg, by="row.names", all.x=TRUE, all.y=TRUE)
colnames(phylo_lunularia_div_classes) <- c("superclass", "phylo_lunularia_div_classes_pos", "phylo_lunularia_div_classes_neg")
phylo_lunularia_div_classes[is.na(phylo_lunularia_div_classes)] <- 0
phylo_lunularia_div_classes$frequency <- apply(X=phylo_lunularia_div_classes[,c(2:ncol(phylo_lunularia_div_classes))], MARGIN=1, FUN=function(x) { sum(x) })
rownames(phylo_lunularia_div_classes) <- phylo_lunularia_div_classes[,1]
phylo_lunularia_div_classes <- phylo_lunularia_div_classes[-which(phylo_lunularia_div_classes$frequency==0),]

pdf(file="plots/phylo_lunularia_ms2_classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(phylo_lunularia_div_classes), phylo_lunularia_div_classes$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(phylo_lunularia_div_classes), "Number of spectra"=phylo_lunularia_div_classes$frequency), file="plots/phylo_lunularia_ms2_classes_sunburst.csv", row.names=FALSE)



# ---------- Histogram ----------
pdf(file="plots/phylo_ms1_phylo_superclass_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(as.matrix(phylo_superclass_list)))
dev.off()



# ---------- Variation partitioning ----------
model_varpart <- varpart(scale(phylo_superclass_list), ~ phylo_mzml_pheno_samples, ~ phylo_mzml_pheno_samples)

# Plot results
pdf(file="plots/phylo_ms1_phylo_superclass_list_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart, Xnames=c("samples","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca <- prcomp(phylo_superclass_list, scale=TRUE, center=TRUE)

# PC-Axes 1+2
pdf(paste("plots/phylo_ms1_phylo_superclass_list_pca12.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 1], model_pca$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC1: ", format(summary(model_pca)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples, bg=phylo_mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 1], model_pca$x[,2], labels=phylo_mzml_names, col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()

# PC-Axes 2+3
pdf(paste("plots/phylo_ms1_phylo_superclass_list_pca23.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca$x[, 2], model_pca$x[,3], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples, bg=phylo_mzml_pheno_colors_samples, cex=2)
grid()
text(model_pca$x[, 2], model_pca$x[,3], labels=phylo_mzml_names, col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
dev.off()



# ---------- Variable Selection ----------
# PLS
phylo_sel_pls <- f.select_features_pls(feat_matrix=phylo_superclass_list, sel_factor=phylo_mzml_pheno_samples, sel_colors=phylo_mzml_pheno_colors, components=3, tune_length=10, quantile_threshold=0.999, plot_roc_filename="plots/phylo_ms1_phylo_superclass_list_select_pls_roc.pdf")
print(paste("Number of selected compounds:", f.count.selected_features(sel_feat=phylo_sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=phylo_superclass_list, sel_feat=phylo_sel_pls$`_selected_variables_`, sample_colors=phylo_mzml_pheno_colors_samples, plot_width=4, plot_height=5, filename="plots/phylo_ms1_phylo_superclass_list_select_pls.pdf", main="")
phylo_sel_pls$`_selected_variables_`
phylo_sel_pls$`_model_r2_`
phylo_sel_pls$`_multiclass_metrics_`



# ############################## Molecular Descriptors ##############################



# ---------- Read SIRIUS identification ----------
# Get SMILES
phylo_sirius_compound_identifications_pos <- read.table(file=paste0("data/phylo_pos_ms2_compound_identifications.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
phylo_sirius_compound_identifications_neg <- read.table(file=paste0("data/phylo_neg_ms2_compound_identifications.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
phylo_sirius_compounds <- c(phylo_sirius_compound_identifications_pos$smiles, phylo_sirius_compound_identifications_neg$smiles)
phylo_sirius_compound_ids <- c(gsub(x=phylo_sirius_compound_identifications_pos$id, pattern=".*_", replacement=""), gsub(x=phylo_sirius_compound_identifications_neg$id, pattern=".*_", replacement=""))

phylo_ms1_def_pos$smiles <- ""
for (i in 1:nrow(phylo_sirius_compound_identifications_pos)) {
	phylo_ms1_def_pos$smiles[which(rownames(phylo_ms1_def_pos)==gsub(x=phylo_sirius_compound_identifications_pos$id[i], pattern=".*_", replacement=""))] <- phylo_sirius_compound_identifications_pos$smiles[i]
}

phylo_ms1_def_neg$smiles <- ""
for (i in 1:nrow(phylo_sirius_compound_identifications_neg)) {
	phylo_ms1_def_neg$smiles[which(rownames(phylo_ms1_def_neg)==gsub(x=phylo_sirius_compound_identifications_neg$id[i], pattern=".*_", replacement=""))] <- phylo_sirius_compound_identifications_neg$smiles[i]
}



# ---------- Calculate molecular phylo_descriptors with RDKit ----------
# Write list of SMILES for rdkit
writeLines(phylo_sirius_compounds, "data/phylo_compound_identifications_smiles.txt")

# Execute rdkit script
system2(command="./rdkit_calculate_descriptors.py", args=c("-i","data/phylo_compound_identifications_smiles.txt", "-o","data/phylo_phylo_descriptors_rdkit.csv"), wait=TRUE)

# Read rdkit_calculate_phylo_descriptors.py output
phylo_rdkit_phylo_descriptors <- read.csv(file="data/phylo_phylo_descriptors_rdkit.csv", header=TRUE, sep=";", fill=TRUE, dec=".", stringsAsFactors=FALSE, check.names=FALSE)
rownames(phylo_rdkit_phylo_descriptors) <- colnames(phylo_sirius_compounds)
phylo_rdkit_phylo_descriptors[,2:ncol(phylo_rdkit_phylo_descriptors)] <- as.data.frame(apply(X=phylo_rdkit_phylo_descriptors[,2:ncol(phylo_rdkit_phylo_descriptors)], MARGIN=2, FUN=function(x) { as.numeric(x) }))
phylo_rdkit_phylo_descriptors <- phylo_rdkit_phylo_descriptors[,2:ncol(phylo_rdkit_phylo_descriptors)]

# Properly format NAs
for (i in 1:nrow(phylo_rdkit_phylo_descriptors)) {
	phylo_rdkit_phylo_descriptors[i, ] <- as.numeric(phylo_rdkit_phylo_descriptors[i, ])
}
for (i in 1:nrow(phylo_rdkit_phylo_descriptors)) {
	phylo_rdkit_phylo_descriptors[i, as.character(phylo_rdkit_phylo_descriptors[i,]) %in% 'list(NULL)'] <- 0
	phylo_rdkit_phylo_descriptors[i, as.character(phylo_rdkit_phylo_descriptors[i,]) %in% 'NA'] <- 0
	phylo_rdkit_phylo_descriptors[i, as.character(phylo_rdkit_phylo_descriptors[i,]) %in% 'NULL'] <- 0
	phylo_rdkit_phylo_descriptors[i, as.character(phylo_rdkit_phylo_descriptors[i,]) %in% 'NaN'] <- 0
	phylo_rdkit_phylo_descriptors[i, as.character(phylo_rdkit_phylo_descriptors[i,]) %in% 'Inf'] <- 0
}
phylo_rdkit_phylo_descriptors <- as.data.frame(phylo_rdkit_phylo_descriptors)



# ---------- Calculate molecular phylo_descriptors from COCONUT library ----------
# Query phylo_descriptors for MetFrag results in positive mode
library(httr)
library(jsonlite)

phylo_coconut_descriptors <- NULL
for (smile in phylo_sirius_compounds) {
	if ( ! is.na(nchar(smile)) ) {
		phylo_coconut_descriptors[[smile]] <- plyr::rbind.fill(phylo_coconut_descriptors[[smile]], tryCatch(f.query_coconut(smile), error=function(x) data.frame(smiles=smile)))
		if (length(phylo_coconut_descriptors[[smile]]) < 1) phylo_coconut_descriptors[[smile]] <- NULL
	}
}

# Calculate general molecular phylo_descriptors from SMILES
number_of_phosphorus <- unlist(lapply(X=lapply(X=phylo_coconut_descriptors, FUN=function(x) { x <- x$molecular_formula }), FUN=function(x) { x <- as.numeric(unlist(Map(function(x) { x <- gsub("(?<=[A-z])(?![0-9])","1",x,perl = TRUE); table(factor(rep(gsub("\\d+", "", x),as.numeric(gsub("\\D+", "", x))), levels="P"))}, regmatches(x, gregexpr("[A-z]+(\\d+)?", x))))) } ))

c_n_ratio <- unlist(lapply(X=phylo_coconut_descriptors, FUN=function(x) { x <- x$number_of_carbons })) / unlist(lapply(X=phylo_coconut_descriptors, FUN=function(x) { x <- x$number_of_nitrogens }))
c_n_ratio[is.infinite(c_n_ratio)] <- 0

c_n_p_ratio <- c_n_ratio / number_of_phosphorus
c_n_p_ratio[is.infinite(c_n_p_ratio)] <- 0
c_n_p_ratio[is.na(c_n_p_ratio)] <- 0

c_p_ratio <- unlist(lapply(X=phylo_coconut_descriptors, FUN=function(x) { x <- x$number_of_carbons })) / number_of_phosphorus
c_p_ratio[is.infinite(c_p_ratio)] <- 0

for (i in 1:length(phylo_coconut_descriptors)) {
	phylo_coconut_descriptors[[i]]$number_of_phosphorus <- number_of_phosphorus[i]
	phylo_coconut_descriptors[[i]]$c_n_ratio <- c_n_ratio[i]
	phylo_coconut_descriptors[[i]]$c_n_p_ratio <- c_n_p_ratio[i]
	phylo_coconut_descriptors[[i]]$c_p_ratio <- c_p_ratio[i]
}

# Build data frame
phylo_descriptors <- NULL
for (i in 1:length(phylo_coconut_descriptors)) {
	phylo_descriptors <- plyr::rbind.fill(phylo_descriptors, as.data.frame(phylo_coconut_descriptors[[i]][1, ! grepl(x=colnames(phylo_coconut_descriptors[[i]]), pattern="(fragments|ertl|absolute_smiles)", perl=TRUE)]))
}

# Filter data frame with only the relevant columns
phylo_coconut_properties <- c("npl_noh_score", "npl_score", "npl_sugar_score", "number_of_carbons",
							  "number_of_nitrogens", "number_of_oxygens", "number_of_rings", "max_number_of_rings", "min_number_of_rings",
							  "total_atom_number", "bond_count",
							  "alogp", "alogp2", "amralogp", "apol", "bpol", "eccentricConnectivityIndexDescriptor",
							  "fmfDescriptor", "fsp3", "fragmentComplexityDescriptor", "hBondAcceptorCount", "hBondDonorCount",
							  "hybridizationRatioDescriptor", "kappaShapeIndex1", "kappaShapeIndex2", "kappaShapeIndex3",
							  "manholdlogp", "petitjeanNumber", "petitjeanShapeTopo", "petitjeanShapeGeom", "lipinskiRuleOf5Failures",
							  "numberSpiroAtoms", "vabcDescriptor", "vertexAdjMagnitude",
							  "zagrebIndex", "tpsaEfficiency", "weinerPolarityNumber",
							  "pfCounts.count",
							  "number_of_phosphorus", "c_n_ratio", "c_n_p_ratio", "c_p_ratio")
phylo_descriptors <- phylo_descriptors[, phylo_coconut_properties]

# Properly format NAs
for (i in 1:nrow(phylo_descriptors)) {
	phylo_descriptors[i, as.character(phylo_descriptors[i,]) %in% 'list(NULL)'] <- 0
	phylo_descriptors[i, as.character(phylo_descriptors[i,]) %in% 'NA'] <- 0
	phylo_descriptors[i, as.character(phylo_descriptors[i,]) %in% 'NULL'] <- 0
	phylo_descriptors[i, as.character(phylo_descriptors[i,]) %in% 'NaN'] <- 0
}

# Remove rows with all NA or all 0
#phylo_descriptors <- phylo_descriptors[rowSums(is.na(phylo_descriptors)) != ncol(phylo_descriptors), ]
#phylo_descriptors <- na.omit(phylo_descriptors)
#phylo_descriptors <- phylo_descriptors[(rowSums(phylo_descriptors) != 0), ]

# Convert to numeric
for (i in 1:nrow(phylo_descriptors)) {
	phylo_descriptors[i, ] <- as.numeric(unlist(phylo_descriptors[i,]))
}



# ---------- Calculate molecular descriptors with CDK ----------
library(rcdk)

phylo_cdk_descriptors <- NULL
for (i in phylo_sirius_compounds) {
	# Get Structure from SMILES
	phylo_cdk_mol = parse.smiles(i)[[1]]
	
	# Get simple measures
	phylo_cdk_atoms = get.atoms(phylo_cdk_mol)
	phylo_cdk_bonds = get.bonds(phylo_cdk_mol)
	
	# Calculate simple measures
	phylo_cdk_num_atoms = as.factor(unlist(lapply(phylo_cdk_atoms, get.symbol)))
	phylo_cdk_num_atoms = tapply(phylo_cdk_num_atoms, phylo_cdk_num_atoms, length)
	numC = as.numeric(phylo_cdk_num_atoms["C"])
	numN = as.numeric(phylo_cdk_num_atoms["N"])
	numP = as.numeric(phylo_cdk_num_atoms["P"])
	numO = as.numeric(phylo_cdk_num_atoms["O"])
	CNRatio = as.numeric(numC / numN)
	
	# Calculate descriptors and restrict to only "constitutional" and "topological"
	phylo_cdk_mol_des_cats = get.desc.categories()
	phylo_cdk_mol_des_names = c(get.desc.names(phylo_cdk_mol_des_cats[3]), get.desc.names(phylo_cdk_mol_des_cats[4]))
	phylo_cdk_mol_des = as.data.frame(eval.desc(phylo_cdk_mol, phylo_cdk_mol_des_names))
	phylo_cdk_descriptors <- plyr::rbind.fill(phylo_cdk_descriptors, cbind(data.frame(numC=numC, numN=numN, numP=numP, numO=numO, CNRatio=CNRatio), phylo_cdk_mol_des))
}

# Properly format NAs and convert to numeric
for (i in 1:nrow(phylo_cdk_descriptors)) {
	phylo_cdk_descriptors[i, as.character(phylo_cdk_descriptors[i,]) %in% 'list(NULL)'] <- 0
	phylo_cdk_descriptors[i, as.character(phylo_cdk_descriptors[i,]) %in% 'NA'] <- 0
	phylo_cdk_descriptors[i, as.character(phylo_cdk_descriptors[i,]) %in% 'NULL'] <- 0
	phylo_cdk_descriptors[i, as.character(phylo_cdk_descriptors[i,]) %in% 'NaN'] <- 0
	phylo_cdk_descriptors[i, ] <- as.numeric(unlist(phylo_cdk_descriptors[i,]))
}



# ---------- Build descriptor list ----------
# Table L: samples x metabolites
phylo_mdes_tab_l <- phylo_comp_list

# Table R: samples x species
phylo_mdes_tab_r <- as.data.frame.matrix(table(rownames(phylo_mdes_tab_l), phylo_mzml_pheno_samples))
rownames(phylo_mdes_tab_r) <- rownames(phylo_mdes_tab_l)

# Table Q: metabolites x traits
phylo_mdes_tab_q <- phylo_cdk_descriptors
phylo_mdes_tab_q[is.na(phylo_mdes_tab_q)] <- 0
phylo_mdes_tab_l <- phylo_comp_list[, which(colnames(phylo_comp_list) %in% c(paste0(rownames(phylo_ms1_def_pos)[which(phylo_ms1_def_pos$smiles != "")], "_pos"), paste0(rownames(phylo_ms1_def_neg)[which(phylo_ms1_def_neg$smiles != "")], "_neg")))]
phylo_mdes_tab_l <- as.data.frame(phylo_mdes_tab_l)

# Perform matrix operation
phylo_mdes_list <- as.data.frame(as.matrix(phylo_mdes_tab_l) %*% as.matrix(phylo_mdes_tab_q))



# ---------- Perform Full dbRDA ----------
attach(phylo_mdes_list)
rm(CNRatio, numC, numN, numO, numP)
phylo_mdes_rda_formula <- formula(paste0("~ 0 + ", paste0(colnames(phylo_mdes_list),collapse=" + ")))
phylo_mdes_rda_y <- data.frame(model.matrix(phylo_mdes_rda_formula))
detach(phylo_mdes_list)

# Calculate overlay of phylo_descriptors on phylo_comp_list distances
phylo_model_mdes_dbrda <- vegan::dbrda(formula=phylo_comp_list ~ ., data=phylo_mdes_rda_y, distance="euclidean")
phylo_model_mdes_dbrda_scores <- vegan::scores(phylo_model_mdes_dbrda, choices=c(1,2))
phylo_model_mdes_dbrda_ef_phylo_descriptors <- envfit(phylo_model_mdes_dbrda, phylo_mdes_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
phylo_model_mdes_dbrda_fit_phylo_descriptors <- data.frame(r2=c(phylo_model_mdes_dbrda_ef_phylo_descriptors$vectors$r,phylo_model_mdes_dbrda_ef_phylo_descriptors$factors$r),
														   pvals=c(phylo_model_mdes_dbrda_ef_phylo_descriptors$vectors$pvals,phylo_model_mdes_dbrda_ef_phylo_descriptors$factors$pvals) )
rownames(phylo_model_mdes_dbrda_fit_phylo_descriptors) <- c(names(phylo_model_mdes_dbrda_ef_phylo_descriptors$vectors$r),names(phylo_model_mdes_dbrda_ef_phylo_descriptors$factors$r))
phylo_model_mdes_dbrda_fit_phylo_descriptors
write.csv(phylo_model_mdes_dbrda_fit_phylo_descriptors, file="plots/phylo_descriptors_dbrda_fit.csv", row.names=TRUE)


# Re-calculate dbRDA with selected variables
phylo_mdes_sel_list <- as.data.frame(as.matrix(phylo_mdes_tab_l) %*% as.matrix(phylo_mdes_tab_q[, rownames(phylo_model_mdes_dbrda_fit_phylo_descriptors)[which(phylo_model_mdes_dbrda_fit_phylo_descriptors$pvals<0.01)]]))
attach(phylo_mdes_sel_list)
phylo_mdes_rda_sel_formula <- formula(paste0("~ 0 + ", paste0(colnames(phylo_mdes_sel_list),collapse=" + ")))
phylo_mdes_rda_sel_y <- data.frame(model.matrix(phylo_mdes_rda_sel_formula))
detach(phylo_mdes_sel_list)
phylo_model_mdes_dbrda_sel <- vegan::dbrda(formula=phylo_comp_list ~ ., data=phylo_mdes_rda_sel_y, distance="euclidean")
phylo_model_mdes_dbrda_sel_scores <- vegan::scores(phylo_model_mdes_dbrda_sel, choices=c(1,2))

# Calculate overlay of phylo_descriptors on metabolite profile distances
phylo_model_mdes_dbrda_sel_ef_phylo_descriptors <- envfit(phylo_model_mdes_dbrda_sel, phylo_mdes_sel_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
phylo_model_mdes_dbrda_sel_fit_phylo_descriptors <- data.frame(r2=c(phylo_model_mdes_dbrda_sel_ef_phylo_descriptors$vectors$r,phylo_model_mdes_dbrda_sel_ef_phylo_descriptors$factors$r),
															   pvals=c(phylo_model_mdes_dbrda_sel_ef_phylo_descriptors$vectors$pvals,phylo_model_mdes_dbrda_sel_ef_phylo_descriptors$factors$pvals) )
rownames(phylo_model_mdes_dbrda_sel_fit_phylo_descriptors) <- c(names(phylo_model_mdes_dbrda_sel_ef_phylo_descriptors$vectors$r),names(phylo_model_mdes_dbrda_sel_ef_phylo_descriptors$factors$r))
phylo_model_mdes_dbrda_sel_fit_phylo_descriptors
write.csv(phylo_model_mdes_dbrda_sel_fit_phylo_descriptors, file="plots/phylo_descriptors_dbrda_sel_fit.csv", row.names=TRUE)


# Plot dbRDA for phylo_descriptors
pdf(file="plots/phylo_descriptors_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(phylo_model_mdes_dbrda_sel_scores$sites[,1])-1, max(phylo_model_mdes_dbrda_sel_scores$sites[,1])+1),
	 ylim=c(min(phylo_model_mdes_dbrda_sel_scores$sites[,2]), max(phylo_model_mdes_dbrda_sel_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(phylo_model_mdes_dbrda_sel)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(phylo_model_mdes_dbrda_sel)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(phylo_model_mdes_dbrda_sel, display="sites", pch=19, col=phylo_mzml_pheno_colors_samples)
text(phylo_model_mdes_dbrda_sel, display="sites", labels=rownames(phylo_comp_list), col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(phylo_model_mdes_dbrda_sel_ef_phylo_descriptors, cex=0.75, p.max=1, col="black")
dev.off()


# Calculate overlay of phylo_comp_list on metabolite profile distances
phylo_model_mdes_dbrda_ef_phylo_comp_list <- envfit(phylo_model_mdes_dbrda, phylo_comp_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
phylo_model_mdes_dbrda_fit_phylo_comp_list <- data.frame(r2=c(phylo_model_mdes_dbrda_ef_phylo_comp_list$vectors$r,phylo_model_mdes_dbrda_ef_phylo_comp_list$factors$r),
														 pvals=c(phylo_model_mdes_dbrda_ef_phylo_comp_list$vectors$pvals,phylo_model_mdes_dbrda_ef_phylo_comp_list$factors$pvals) )
rownames(phylo_model_mdes_dbrda_fit_phylo_comp_list) <- c(names(phylo_model_mdes_dbrda_ef_phylo_comp_list$vectors$r),names(phylo_model_mdes_dbrda_ef_phylo_comp_list$factors$r))
phylo_model_mdes_dbrda_fit_phylo_comp_list
write.csv(phylo_model_mdes_dbrda_fit_phylo_comp_list, file="plots/phylo_comp_list_dbrda_fit.csv", row.names=TRUE)

# Plot dbRDA for phylo_descriptors
pdf(file="plots/phylo_comp_list_dbrda_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(phylo_model_mdes_dbrda_scores$sites[,1])-1, max(phylo_model_mdes_dbrda_scores$sites[,1])+1),
	 ylim=c(min(phylo_model_mdes_dbrda_scores$sites[,2]), max(phylo_model_mdes_dbrda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(phylo_model_mdes_dbrda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(phylo_model_mdes_dbrda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(phylo_model_mdes_dbrda, display="sites", pch=19, col=phylo_mzml_pheno_colors_samples)
text(phylo_model_mdes_dbrda, display="sites", labels=rownames(phylo_comp_list), col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(phylo_model_mdes_dbrda_ef_phylo_comp_list, cex=0.75, p.max=1, col="black")
dev.off()


# Calculate overlay of phylo_class_list on metabolite profile distances
phylo_model_mdes_dbrda_ef_phylo_class_list <- envfit(phylo_model_mdes_dbrda, phylo_class_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
phylo_model_mdes_dbrda_fit_phylo_class_list <- data.frame(r2=c(phylo_model_mdes_dbrda_ef_phylo_class_list$vectors$r,phylo_model_mdes_dbrda_ef_phylo_class_list$factors$r),
														  pvals=c(phylo_model_mdes_dbrda_ef_phylo_class_list$vectors$pvals,phylo_model_mdes_dbrda_ef_phylo_class_list$factors$pvals) )
rownames(phylo_model_mdes_dbrda_fit_phylo_class_list) <- c(names(phylo_model_mdes_dbrda_ef_phylo_class_list$vectors$r),names(phylo_model_mdes_dbrda_ef_phylo_class_list$factors$r))
phylo_model_mdes_dbrda_fit_phylo_class_list
write.csv(phylo_model_mdes_dbrda_fit_phylo_class_list, file="plots/phylo_class_list_dbrda_fit.csv", row.names=TRUE)

# Plot dbRDA for phylo_descriptors
pdf(file="plots/phylo_class_list_dbrda_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(phylo_model_mdes_dbrda_scores$sites[,1])-1, max(phylo_model_mdes_dbrda_scores$sites[,1])+1),
	 ylim=c(min(phylo_model_mdes_dbrda_scores$sites[,2]), max(phylo_model_mdes_dbrda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(phylo_model_mdes_dbrda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(phylo_model_mdes_dbrda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(phylo_model_mdes_dbrda, display="sites", pch=19, col=phylo_mzml_pheno_colors_samples)
text(phylo_model_mdes_dbrda, display="sites", labels=rownames(phylo_comp_list), col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(phylo_model_mdes_dbrda_ef_phylo_class_list, cex=0.75, p.max=1, col="black")
dev.off()


# Calculate overlay of phylo_superclass_list on metabolite profile distances
phylo_model_mdes_dbrda_ef_phylo_superclass_list <- envfit(phylo_model_mdes_dbrda, phylo_superclass_list, perm=10000)

# Goodness of fit statistic: Squared correlation coefficient
phylo_model_mdes_dbrda_fit_phylo_superclass_list <- data.frame(r2=c(phylo_model_mdes_dbrda_ef_phylo_superclass_list$vectors$r,phylo_model_mdes_dbrda_ef_phylo_superclass_list$factors$r),
															   pvals=c(phylo_model_mdes_dbrda_ef_phylo_superclass_list$vectors$pvals,phylo_model_mdes_dbrda_ef_phylo_superclass_list$factors$pvals) )
rownames(phylo_model_mdes_dbrda_fit_phylo_superclass_list) <- c(names(phylo_model_mdes_dbrda_ef_phylo_superclass_list$vectors$r),names(phylo_model_mdes_dbrda_ef_phylo_superclass_list$factors$r))
phylo_model_mdes_dbrda_fit_phylo_superclass_list
write.csv(phylo_model_mdes_dbrda_fit_phylo_superclass_list, file="plots/phylo_superclass_list_dbrda_fit.csv", row.names=TRUE)

# Plot dbRDA for phylo_descriptors
pdf(file="plots/phylo_superclass_list_dbrda_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(phylo_model_mdes_dbrda_scores$sites[,1])-1, max(phylo_model_mdes_dbrda_scores$sites[,1])+1),
	 ylim=c(min(phylo_model_mdes_dbrda_scores$sites[,2]), max(phylo_model_mdes_dbrda_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(phylo_model_mdes_dbrda)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(phylo_model_mdes_dbrda)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(phylo_model_mdes_dbrda, display="sites", pch=19, col=phylo_mzml_pheno_colors_samples)
text(phylo_model_mdes_dbrda, display="sites", labels=rownames(phylo_comp_list), col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(phylo_model_mdes_dbrda_ef_phylo_superclass_list, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- dbRDA of phylo_mdes_list on phylo_comp_list with stepwise forward selection ----------
attach(phylo_mdes_list)
phylo_mdes_rda_comp_formula <- formula(paste0("~ 0 + ", paste0(colnames(phylo_mdes_list),collapse=" + ")))
phylo_mdes_rda_comp_y <- data.frame(model.matrix(phylo_mdes_rda_comp_formula))
detach(phylo_mdes_list)

# Model with intercept only
phylo_model_0_mdes_rda_comp <- vegan::dbrda(formula=phylo_comp_list ~ 1, data=phylo_mdes_rda_comp_y, distance="euclidean")

# Model with all explanatory variables
phylo_model_1_mdes_rda_comp <- vegan::dbrda(formula=phylo_comp_list ~ ., data=phylo_mdes_rda_comp_y, distance="euclidean")

# Stepwise forward selection
phylo_model_step_mdes_rda_comp <- ordistep(object=phylo_model_0_mdes_rda_comp, scope=formula(phylo_model_1_mdes_rda_comp), Pin=0.1, Pout=0.5, direction="forward", perm.max=100000)
phylo_model_step_mdes_rda_comp_scores <- vegan::scores(phylo_model_step_mdes_rda_comp)

# dbRDA with selected model by permutation tests in constrained ordination
phylo_model_mdes_rda_comp <- vegan::dbrda(formula=as.formula(phylo_model_step_mdes_rda_comp$terms), data=phylo_mdes_rda_comp_y, distance="euclidean")
phylo_model_mdes_rda_comp_ef_formula <- update(as.formula(phylo_model_step_mdes_rda_comp$terms), phylo_model_mdes_rda_comp ~ .)
phylo_model_mdes_rda_comp_ef_factors <- as.factor(sapply(strsplit(as.character(phylo_model_mdes_rda_comp_ef_formula)[[3]], "\\+"), function(x) { x <- gsub("(\\`|^ | $)","",x) }))
phylo_model_mdes_rda_comp_ef <- envfit(formula=phylo_model_mdes_rda_comp_ef_formula, data=phylo_mdes_rda_comp_y, perm=10000)
phylo_model_mdes_rda_comp_scores <- vegan::scores(phylo_model_mdes_rda_comp)

# Goodness of fit statistic: Squared correlation coefficient
phylo_model_mdes_rda_comp_fit <- data.frame(r2=c(phylo_model_mdes_rda_comp_ef$vectors$r,phylo_model_mdes_rda_comp_ef$factors$r),
											pvals=c(phylo_model_mdes_rda_comp_ef$vectors$pvals,phylo_model_mdes_rda_comp_ef$factors$pvals) )
rownames(phylo_model_mdes_rda_comp_fit) <- c(names(phylo_model_mdes_rda_comp_ef$vectors$r),names(phylo_model_mdes_rda_comp_ef$factors$r))
phylo_model_mdes_rda_comp_fit
write.csv(phylo_model_mdes_rda_comp_fit, file="plots/phylo_ms1_comp_list_dbrda_sel_fit.csv", row.names=TRUE)

# Plot results
pdf(file="plots/phylo_ms1_comp_list_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(phylo_model_mdes_rda_comp_scores$sites[,1])-1, max(phylo_model_mdes_rda_comp_scores$sites[,1])+1),
	 ylim=c(min(phylo_model_mdes_rda_comp_scores$sites[,2]), max(phylo_model_mdes_rda_comp_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(phylo_model_mdes_rda_comp)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(phylo_model_mdes_rda_comp)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(phylo_model_mdes_rda_comp, display="sites", pch=19, col=phylo_mzml_pheno_colors_samples)
text(phylo_model_mdes_rda_comp, display="sites", labels=rownames(phylo_comp_list), col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(phylo_model_mdes_rda_comp_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- dbRDA of phylo_mdes_list on phylo_class_list with stepwise forward selection ----------
attach(phylo_mdes_list)
phylo_mdes_rda_class_formula <- formula(paste0("~ 0 + ", paste0(colnames(phylo_mdes_list),collapse=" + ")))
phylo_mdes_rda_class_y <- data.frame(model.matrix(phylo_mdes_rda_class_formula))
detach(phylo_mdes_list)

# Model with intercept only
phylo_model_0_mdes_rda_class <- vegan::dbrda(formula=phylo_class_list ~ 1, data=phylo_mdes_rda_class_y, distance="euclidean")

# Model with all explanatory variables
phylo_model_1_mdes_rda_class <- vegan::dbrda(formula=phylo_class_list ~ ., data=phylo_mdes_rda_class_y, distance="euclidean")

# Stepwise forward selection
phylo_model_step_mdes_rda_class <- ordistep(object=phylo_model_0_mdes_rda_class, scope=formula(phylo_model_1_mdes_rda_class), Pin=0.1, Pout=0.5, direction="forward", perm.max=100000)
phylo_model_step_mdes_rda_class_scores <- vegan::scores(phylo_model_step_mdes_rda_class)

# dbRDA with selected model by permutation tests in constrained ordination
phylo_model_mdes_rda_class <- vegan::dbrda(formula=as.formula(phylo_model_step_mdes_rda_class$terms), data=phylo_mdes_rda_class_y, distance="euclidean")
phylo_model_mdes_rda_class_ef_formula <- update(as.formula(phylo_model_step_mdes_rda_class$terms), phylo_model_mdes_rda_class ~ .)
phylo_model_mdes_rda_class_ef_factors <- as.factor(sapply(strsplit(as.character(phylo_model_mdes_rda_class_ef_formula)[[3]], "\\+"), function(x) { x <- gsub("(\\`|^ | $)","",x) }))
phylo_model_mdes_rda_class_ef <- envfit(formula=phylo_model_mdes_rda_class_ef_formula, data=phylo_mdes_rda_class_y, perm=10000)
phylo_model_mdes_rda_class_scores <- vegan::scores(phylo_model_mdes_rda_class)

# Goodness of fit statistic: Squared correlation coefficient
phylo_model_mdes_rda_class_fit <- data.frame(r2=c(phylo_model_mdes_rda_class_ef$vectors$r,phylo_model_mdes_rda_class_ef$factors$r),
											 pvals=c(phylo_model_mdes_rda_class_ef$vectors$pvals,phylo_model_mdes_rda_class_ef$factors$pvals) )
rownames(phylo_model_mdes_rda_class_fit) <- c(names(phylo_model_mdes_rda_class_ef$vectors$r),names(phylo_model_mdes_rda_class_ef$factors$r))
phylo_model_mdes_rda_class_fit
write.csv(phylo_model_mdes_rda_class_fit, file="plots/phylo_ms2_class_list_dbrda_sel_fit.csv", row.names=TRUE)

# Plot results
pdf(file="plots/phylo_ms2_class_list_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(phylo_model_mdes_rda_class_scores$sites[,1])-1, max(phylo_model_mdes_rda_class_scores$sites[,1])+1),
	 ylim=c(min(phylo_model_mdes_rda_class_scores$sites[,2]), max(phylo_model_mdes_rda_class_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(phylo_model_mdes_rda_class)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(phylo_model_mdes_rda_class)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(phylo_model_mdes_rda_class, display="sites", pch=19, col=phylo_mzml_pheno_colors_samples)
text(phylo_model_mdes_rda_class, display="sites", labels=rownames(phylo_class_list), col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(phylo_model_mdes_rda_class_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- dbRDA of phylo_mdes_list on phylo_superclass_list with stepwise forward selection ----------
attach(phylo_mdes_list)
phylo_mdes_rda_superclass_formula <- formula(paste0("~ 0 + ", paste0(colnames(phylo_mdes_list),collapse=" + ")))
phylo_mdes_rda_superclass_y <- data.frame(model.matrix(phylo_mdes_rda_superclass_formula))
detach(phylo_mdes_list)

# Model with intercept only
phylo_model_0_mdes_rda_superclass <- vegan::dbrda(formula=phylo_superclass_list ~ 1, data=phylo_mdes_rda_superclass_y, distance="euclidean")

# Model with all explanatory variables
phylo_model_1_mdes_rda_superclass <- vegan::dbrda(formula=phylo_superclass_list ~ ., data=phylo_mdes_rda_superclass_y, distance="euclidean")

# Stepwise forward selection
phylo_model_step_mdes_rda_superclass <- ordistep(object=phylo_model_0_mdes_rda_superclass, scope=formula(phylo_model_1_mdes_rda_superclass), Pin=0.1, Pout=0.5, direction="forward", perm.max=100000)
phylo_model_step_mdes_rda_superclass_scores <- vegan::scores(phylo_model_step_mdes_rda_superclass)

# dbRDA with selected model by permutation tests in constrained ordination
phylo_model_mdes_rda_superclass <- vegan::dbrda(formula=as.formula(phylo_model_step_mdes_rda_superclass$terms), data=phylo_mdes_rda_superclass_y, distance="euclidean")
phylo_model_mdes_rda_superclass_ef_formula <- update(as.formula(phylo_model_step_mdes_rda_superclass$terms), phylo_model_mdes_rda_superclass ~ .)
phylo_model_mdes_rda_superclass_ef_factors <- as.factor(sapply(strsplit(as.character(phylo_model_mdes_rda_superclass_ef_formula)[[3]], "\\+"), function(x) { x <- gsub("(\\`|^ | $)","",x) }))
phylo_model_mdes_rda_superclass_ef <- envfit(formula=phylo_model_mdes_rda_superclass_ef_formula, data=phylo_mdes_rda_superclass_y, perm=10000)
phylo_model_mdes_rda_superclass_scores <- vegan::scores(phylo_model_mdes_rda_superclass)

# Goodness of fit statistic: Squared correlation coefficient
phylo_model_mdes_rda_superclass_fit <- data.frame(r2=c(phylo_model_mdes_rda_superclass_ef$vectors$r,phylo_model_mdes_rda_superclass_ef$factors$r),
												  pvals=c(phylo_model_mdes_rda_superclass_ef$vectors$pvals,phylo_model_mdes_rda_superclass_ef$factors$pvals) )
rownames(phylo_model_mdes_rda_superclass_fit) <- c(names(phylo_model_mdes_rda_superclass_ef$vectors$r),names(phylo_model_mdes_rda_superclass_ef$factors$r))
phylo_model_mdes_rda_superclass_fit
write.csv(phylo_model_mdes_rda_superclass_fit, file="plots/phylo_ms2_superclass_list_dbrda_sel_fit.csv", row.names=TRUE)

# Plot results
pdf(file="plots/phylo_ms2_superclass_list_dbrda_sel_fit.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
plot(0, 0, xlim=c(min(phylo_model_mdes_rda_superclass_scores$sites[,1])-1, max(phylo_model_mdes_rda_superclass_scores$sites[,1])+1),
	 ylim=c(min(phylo_model_mdes_rda_superclass_scores$sites[,2]), max(phylo_model_mdes_rda_superclass_scores$sites[,2])),
	 xlab=paste0("dbRDA1 (",round(as.data.frame(summary(phylo_model_mdes_rda_superclass)$cont$importance)["Proportion Explained",c(1)]*100,2),"%)"),
	 ylab=paste0("dbRDA2 (",round(as.data.frame(summary(phylo_model_mdes_rda_superclass)$cont$importance)["Proportion Explained",c(2)]*100,2),"%)"),
	 main="dbRDA")
points(phylo_model_mdes_rda_superclass, display="sites", pch=19, col=phylo_mzml_pheno_colors_samples)
text(phylo_model_mdes_rda_superclass, display="sites", labels=rownames(phylo_superclass_list), col=phylo_mzml_pheno_colors_samples, pos=3, cex=0.5)
plot(phylo_model_mdes_rda_superclass_ef, cex=0.75, p.max=1, col="black")
dev.off()



# ---------- PLS ----------
# PLS
sel_pls <- f.select_features_pls(feat_matrix=phylo_mdes_list, sel_factor=phylo_mzml_pheno_samples, sel_colors=mzml_pheno_colors, components=2, tune_length=10, quantile_threshold=0.995, plot_roc_filename="plots/phylo_descriptors_select_pls_roc.pdf")
print(paste("Number of selected phylo_descriptors:", f.count.selected_features(sel_feat=sel_pls$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=phylo_mdes_list, sel_feat=sel_pls$`_selected_variables_`, sample_colors=phylo_mzml_pheno_colors_samples, plot_width=4, plot_height=5, filename="plots/phylo_descriptors_select_pls.pdf", main="")
sel_pls$`_multiclass_metrics_`
sel_pls$`_model_r2_`



# ############################## Chemotaxonomy Outgroup ##############################



# ---------- Taxonomic tree from genetic data ----------
# Read phylogenetic tree
phylo_outgroup_tree <- read.tree("data/trnLF_outgroup_newick.txt")
phylo_outgroup_tree <- chronos(phylo_outgroup_tree)
phylo_outgroup_heights <- phylo_outgroup_tree$edge.length
phylo_outgroup_dend <- as.dendrogram(phylo_outgroup_tree)
labels(phylo_outgroup_dend) <- c("R.sorocarpa", "R.glauca", "R.warnstorfii", "L. cruciata")

# Distance matrix of phylogenetic tree using Bray-Curtis
phylo_outgroup_dist <- cophenetic.phylo(phylo_outgroup_tree)
phylo_outgroup_dist <- vegdist(phylo_outgroup_dist, method="bray")

# Hierarchical clustering
phylo_outgroup_hclust <- hclust(phylo_outgroup_dist, method="complete")



# ---------- Chemotaxonomic tree for feat_list ----------
# Merge feat_list for species from samples
phylo_med_feat_list <- NULL
for (i in unique(phylo_mzml_pheno_samples_pos)) phylo_med_feat_list <- rbind(phylo_med_feat_list, apply(X=phylo_feat_list[which(phylo_mzml_pheno_samples_pos==i),], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_med_feat_list) <- unique(phylo_mzml_pheno_samples_pos)

# Distance matrix of feat_list
phylo_feat_dist <- vegdist(phylo_med_feat_list, method="euclidean")

# Hierarchical clustering
phylo_feat_hclust <- hclust(phylo_feat_dist, method="complete")

# Optimal order
#phylo_feat_opti <- order.optimal(phylo_feat_dist, phylo_feat_hclust$merge)
phylo_feat_oclust <- phylo_feat_hclust
#phylo_feat_oclust$merge <- phylo_feat_opti$merge
#phylo_feat_oclust$order <- phylo_feat_opti$order

# Dendrogram
phylo_feat_dend <- as.dendrogram(phylo_feat_oclust)



# ---------- Chemotaxonomic tree for compound list ----------
# Merge comp_list for species from samples
phylo_med_comp_list <- NULL
for (i in unique(phylo_mzml_pheno_samples_pos)) phylo_med_comp_list <- rbind(phylo_med_comp_list, apply(X=phylo_comp_list[phylo_mzml_pheno_samples_pos==i,], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_med_comp_list) <- unique(phylo_mzml_pheno_samples_pos)

# Distance matrix of feat_list using Bray-Curtis
phylo_comp_dist <- vegdist(phylo_med_comp_list, method="euclidean")

# Hierarchical clustering
phylo_comp_hclust <- hclust(phylo_comp_dist, method="complete")

# Optimal order
#phylo_comp_opti <- order.optimal(phylo_comp_dist, phylo_comp_hclust$merge)
phylo_comp_oclust <- phylo_comp_hclust
#phylo_comp_oclust$merge <- phylo_comp_opti$merge
#phylo_comp_oclust$order <- phylo_comp_opti$order

# Dendrogram
phylo_comp_dend <- as.dendrogram(phylo_comp_oclust)



# ---------- Chemotaxonomic tree for class list ----------
# Merge class_list for species from samples
phylo_med_class_list <- NULL
for (i in unique(phylo_mzml_pheno_samples_pos)) phylo_med_class_list <- rbind(phylo_med_class_list, apply(X=phylo_class_list[phylo_mzml_pheno_samples_pos==i,], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_med_class_list) <- unique(phylo_mzml_pheno_samples_pos)

# Distance matrix of feat_list using Bray-Curtis
phylo_class_dist <- vegdist(phylo_med_class_list, method="euclidean")

# Hierarchical clustering
phylo_class_hclust <- hclust(phylo_class_dist, method="complete")

# Optimal order
#phylo_class_opti <- order.optimal(phylo_class_dist, phylo_class_hclust$merge)
phylo_class_oclust <- phylo_class_hclust
#phylo_class_oclust$merge <- phylo_class_opti$merge
#phylo_class_oclust$order <- phylo_class_opti$order

# Dendrogram
phylo_class_dend <- as.dendrogram(phylo_class_oclust)



# ---------- Chemotaxonomic tree for superclass list ----------
# Merge superclass_list for species from samples
phylo_med_superclass_list <- NULL
for (i in unique(phylo_mzml_pheno_samples_pos)) phylo_med_superclass_list <- rbind(phylo_med_superclass_list, apply(X=phylo_superclass_list[phylo_mzml_pheno_samples_pos==i,], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_med_superclass_list) <- unique(phylo_mzml_pheno_samples_pos)

# Distance matrix of feat_list using Bray-Curtis
phylo_superclass_dist <- vegdist(phylo_med_superclass_list, method="euclidean")

# Hierarchical clustering
phylo_superclass_hclust <- hclust(phylo_superclass_dist, method="complete")

# Optimal order
#phylo_superclass_opti <- order.optimal(phylo_superclass_dist, phylo_superclass_hclust$merge)
phylo_superclass_oclust <- phylo_superclass_hclust
#phylo_superclass_oclust$merge <- phylo_superclass_opti$merge
#phylo_superclass_oclust$order <- phylo_superclass_opti$order

# Dendrogram
phylo_superclass_dend <- as.dendrogram(phylo_superclass_oclust)



# ---------- Chemotaxonomic tree for mdes list ----------
# Merge mdes_list for species from samples
phylo_med_mdes_list <- NULL
for (i in unique(phylo_mzml_pheno_samples_pos)) phylo_med_mdes_list <- rbind(phylo_med_mdes_list, apply(X=phylo_mdes_list[phylo_mzml_pheno_samples_pos==i,], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(phylo_med_mdes_list) <- unique(phylo_mzml_pheno_samples_pos)

# Distance matrix of feat_list using Bray-Curtis
phylo_mdes_dist <- vegdist(phylo_med_mdes_list, method="bray")

# Hierarchical clustering
phylo_mdes_hclust <- hclust(phylo_mdes_dist, method="complete")

# Optimal order
#phylo_mdes_opti <- order.optimal(phylo_mdes_dist, phylo_mdes_hclust$merge)
phylo_mdes_oclust <- phylo_mdes_hclust
#phylo_mdes_oclust$merge <- phylo_mdes_opti$merge
#phylo_mdes_oclust$order <- phylo_mdes_opti$order

# Dendrogram
phylo_mdes_dend <- as.dendrogram(phylo_mdes_oclust)



# ---------- Plot chemotaxonomic trees ----------
# Plot chemotaxonomic trees
pdf("plots/chemotax_phylo_outgroup_trees.pdf", encoding="ISOLatin1", pointsize=8, width=5*2, height=2.8, family="Helvetica")
par(mfrow=c(1,5), mar=c(1,1,2,1), cex=1.0)
labels(phylo_outgroup_hclust) <- c("Lunularia cruciata", "R. sorocarpa", "R. glauca", "R. warnstorfii")
phylo_max <- max(phylo_outgroup_dist)*1.2; plot(as.phylo(phylo_outgroup_hclust), type="phylogram", direction="leftwards", x.lim=c(0,phylo_max), label.offset=phylo_max/30, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=c("seagreen3","seagreen3","seagreen3","forestgreen"), font=4, main="")
edgelabels(text=round(node.depth.edgelength(as.phylo(phylo_outgroup_hclust))[-5],2), col="black", bg="white", frame="none", adj=c(0.5,-0.5), cex=0.8)
mtext(text="(a)", adj=0, line=0.5, font=2, cex=1.2)
phylo_max <- max(as.dist(phylo_comp_oclust))*1.2; plot(as.phylo(phylo_comp_oclust), type="phylogram", direction="rightwards", x.lim=c(0,phylo_max), label.offset=phylo_max/30, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=c("forestgreen","seagreen3","seagreen3","seagreen3"), font=4, main="")
edgelabels(text=round(node.depth.edgelength(as.phylo(phylo_comp_oclust))[-5],2), col="black", bg="white", frame="none", adj=c(0.5,-0.5), cex=0.8)
mtext(text="(b)", adj=0, line=0.5, font=2, cex=1.2)
phylo_max <- max(as.dist(phylo_class_oclust))*1.2; plot(as.phylo(phylo_class_oclust), type="phylogram", direction="rightwards", x.lim=c(0,phylo_max), label.offset=phylo_max/30, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=c("forestgreen","seagreen3","seagreen3","seagreen3"), font=4, main="")
edgelabels(text=round(node.depth.edgelength(as.phylo(phylo_class_oclust))[-5],2), col="black", bg="white", frame="none", adj=c(0.5,-0.5), cex=0.8)
mtext(text="(c)", adj=0, line=0.5, font=2, cex=1.2)
phylo_max <- max(as.dist(phylo_superclass_oclust))*1.2; plot(as.phylo(phylo_superclass_oclust), type="phylogram", direction="rightwards", x.lim=c(0,phylo_max), label.offset=phylo_max/30, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=c("forestgreen","seagreen3","seagreen3","seagreen3"), font=4, main="")
edgelabels(text=round(node.depth.edgelength(as.phylo(phylo_superclass_oclust))[-5],2), col="black", bg="white", frame="none", adj=c(0.5,-0.5), cex=0.8)
mtext(text="(d)", adj=0, line=0.5, font=2, cex=1.2)
phylo_max <- max(as.dist(phylo_mdes_oclust))*1.2; plot(as.phylo(phylo_mdes_oclust), type="phylogram", direction="rightwards", x.lim=c(0,phylo_max), label.offset=phylo_max/30, use.edge.length=TRUE, show.tip.label=TRUE, tip.color=c("forestgreen","seagreen3","seagreen3","seagreen3"), font=4, main="")
edgelabels(text=round(node.depth.edgelength(as.phylo(phylo_mdes_oclust))[-5],2), col="black", bg="white", frame="none", adj=c(0.5,-0.5), cex=0.8)
mtext(text="(e)", adj=0, line=0.5, font=2, cex=1.2)
dev.off()



# ---------- Clustering metrics ----------
phylo_model_tree_metrics <- NULL
phylo_model_tree_metrics <- rbind(phylo_model_tree_metrics, c("mantel"=round(mantel(xdis=phylo_outgroup_dist, ydis=phylo_comp_dist, method="pearson", permutations=10000)$statistic, 3),
															  "cor"=cor(phylo_outgroup_dist, phylo_comp_dist, method="pearson"),
															  "cop"=cor_cophenetic(hclust(phylo_outgroup_dist), hclust(phylo_comp_dist), method="pearson"),
															  "robinson-foulds"=RF.dist(as.phylo(phylo_outgroup_dend), as.phylo(phylo_comp_oclust), normalize=TRUE, check.labels=FALSE, rooted=TRUE)
))
phylo_model_tree_metrics <- rbind(phylo_model_tree_metrics, c("mantel"=round(mantel(xdis=phylo_outgroup_dist, ydis=phylo_class_dist, method="pearson", permutations=10000)$statistic, 3),
															  "cor"=cor(phylo_outgroup_dist, phylo_class_dist, method="pearson"),
															  "cop"=cor_cophenetic(hclust(phylo_outgroup_dist), hclust(phylo_class_dist), method="pearson"),
															  "robinson-foulds"=RF.dist(as.phylo(phylo_outgroup_dend), as.phylo(phylo_class_oclust), normalize=TRUE, check.labels=FALSE, rooted=TRUE)
))
phylo_model_tree_metrics <- rbind(phylo_model_tree_metrics, c("mantel"=round(mantel(xdis=phylo_outgroup_dist, ydis=phylo_superclass_dist, method="pearson", permutations=10000)$statistic, 3),
															  "cor"=cor(phylo_outgroup_dist, phylo_superclass_dist, method="pearson"),
															  "cop"=cor_cophenetic(hclust(phylo_outgroup_dist), hclust(phylo_superclass_dist), method="pearson"),
															  "robinson-foulds"=RF.dist(as.phylo(phylo_outgroup_dend), as.phylo(phylo_superclass_oclust), normalize=TRUE, check.labels=FALSE, rooted=TRUE)
))
phylo_model_tree_metrics <- rbind(phylo_model_tree_metrics, c("mantel"=round(mantel(xdis=phylo_outgroup_dist, ydis=phylo_mdes_dist, method="pearson", permutations=10000)$statistic, 3),
															  "cor"=cor(phylo_outgroup_dist, phylo_mdes_dist, method="pearson"),
															  "cop"=cor_cophenetic(hclust(phylo_outgroup_dist), hclust(phylo_mdes_dist), method="pearson"),
															  "robinson-foulds"=RF.dist(as.phylo(phylo_outgroup_dend), as.phylo(phylo_mdes_oclust), normalize=TRUE, check.labels=FALSE, rooted=TRUE)
))
rownames(phylo_model_tree_metrics) <- c("comp", "class", "superclass", "mdes")
write.csv(phylo_model_tree_metrics, file="plots/chemotax_phylo_outgroup_trees.csv", row.names=TRUE)



# ############################## MetaboLights Export ##############################



# ---------- Ingroup ----------
riccia_ingroup_table <- riccia_id
riccia_ingroup_xcms <- gsub(x=riccia_ingroup_table$XCMS.name, pattern="_.*", replacement="")
riccia_ingroup_mode <- gsub(x=riccia_ingroup_table$XCMS.name, pattern=".*_", replacement="")
riccia_ingroup_smiles <- riccia_ingroup_table$SMILES
riccia_ingroup_name <- riccia_ingroup_table$Compound.Name
for (i in 1:length(riccia_ingroup_smiles)) {
	if ((riccia_ingroup_smiles[i] == "-") | (riccia_ingroup_smiles[i] == "n.d")) {
		riccia_ingroup_smiles[i] <- NA
	}
	if ((riccia_ingroup_name[i] == "-") | (riccia_ingroup_name[i] == "n.d")) {
		riccia_ingroup_name[i] <- NA
	}
	if (is.na(riccia_ingroup_smiles[i])) {
		riccia_ingroup_xcms[i] <- NA
		riccia_ingroup_mode[i] <- NA
	}
}

riccia_ingroup_xcms <- na.omit(riccia_ingroup_xcms)
riccia_ingroup_mode <- na.omit(riccia_ingroup_mode)
riccia_ingroup_smiles <- na.omit(riccia_ingroup_smiles)
riccia_ingroup_name <- na.omit(riccia_ingroup_name)

# Annotate feat_list in neg-mode
riccia_ingroup_maf_neg <- ms1_def_neg
f.export_maf(riccia_ingroup_maf_neg, "data/ingroup_maf_neg.tsv")
f.annotate_maf_classes(maf_input="data/ingroup_maf_neg.tsv", maf_output="data/ingroup_maf_neg_annotated_classes.tsv")
f.annotate_maf_compounds(maf_input="data/ingroup_maf_neg_annotated_classes.tsv", maf_output="data/ingroup_maf_neg_annotated_compounds.tsv", polarity="neg", xcms_id=riccia_ingroup_xcms, pol_mode=riccia_ingroup_mode, smiles=riccia_ingroup_smiles, names=riccia_ingroup_name)

# Annotate feat_list in pos-mode
riccia_ingroup_maf_pos <- ms1_def_pos
f.export_maf(riccia_ingroup_maf_pos, "data/ingroup_maf_pos.tsv")
f.annotate_maf_classes(maf_input="data/ingroup_maf_pos.tsv", maf_output="data/ingroup_maf_pos_annotated_classes.tsv")
f.annotate_maf_compounds(maf_input="data/ingroup_maf_pos_annotated_classes.tsv", maf_output="data/ingroup_maf_pos_annotated_compounds.tsv", polarity="pos", xcms_id=riccia_ingroup_xcms, pol_mode=riccia_ingroup_mode, smiles=riccia_ingroup_smiles, names=riccia_ingroup_name)



# ---------- Outgroup ----------
riccia_outgroup_table <- riccia_id
riccia_outgroup_xcms <- gsub(x=riccia_outgroup_table$XCMS.name, pattern="_.*", replacement="")
riccia_outgroup_mode <- gsub(x=riccia_outgroup_table$XCMS.name, pattern=".*_", replacement="")
riccia_outgroup_smiles <- riccia_outgroup_table$SMILES
riccia_outgroup_name <- riccia_outgroup_table$Compound.Name
for (i in 1:length(riccia_outgroup_smiles)) {
	if ((riccia_outgroup_smiles[i] == "-") | (riccia_outgroup_smiles[i] == "n.d")) {
		riccia_outgroup_smiles[i] <- NA
	}
	if ((riccia_outgroup_name[i] == "-") | (riccia_outgroup_name[i] == "n.d")) {
		riccia_outgroup_name[i] <- NA
	}
	if (is.na(riccia_outgroup_smiles[i])) {
		riccia_outgroup_xcms[i] <- NA
		riccia_outgroup_mode[i] <- NA
	}
}

riccia_outgroup_xcms <- na.omit(riccia_outgroup_xcms)
riccia_outgroup_mode <- na.omit(riccia_outgroup_mode)
riccia_outgroup_smiles <- na.omit(riccia_outgroup_smiles)
riccia_outgroup_name <- na.omit(riccia_outgroup_name)

# Annotate feat_list in neg-mode
riccia_outgroup_maf_neg <- ms1_def_neg
f.export_maf(riccia_outgroup_maf_neg, "data/outgroup_maf_neg.tsv")
f.annotate_maf_classes(maf_input="data/outgroup_maf_neg.tsv", maf_output="data/outgroup_maf_neg_annotated_classes.tsv")
f.annotate_maf_compounds(maf_input="data/outgroup_maf_neg_annotated_classes.tsv", maf_output="data/outgroup_maf_neg_annotated_compounds.tsv", polarity="neg", xcms_id=riccia_outgroup_xcms, pol_mode=riccia_outgroup_mode, smiles=riccia_outgroup_smiles, names=riccia_outgroup_name)

# Annotate feat_list in pos-mode
riccia_outgroup_maf_pos <- ms1_def_pos
f.export_maf(riccia_outgroup_maf_pos, "data/outgroup_maf_pos.tsv")
f.annotate_maf_classes(maf_input="data/outgroup_maf_pos.tsv", maf_output="data/outgroup_maf_pos_annotated_classes.tsv")
f.annotate_maf_compounds(maf_input="data/outgroup_maf_pos_annotated_classes.tsv", maf_output="data/outgroup_maf_pos_annotated_compounds.tsv", polarity="pos", xcms_id=riccia_outgroup_xcms, pol_mode=riccia_outgroup_mode, smiles=riccia_outgroup_smiles, names=riccia_outgroup_name)


