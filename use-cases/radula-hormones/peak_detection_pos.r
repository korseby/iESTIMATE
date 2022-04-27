


# ############################## MS1 ##############################



# ---------- MS1 Preparations ----------
# MS1 variables
polarity <- "positive"
pol <- substr(x=polarity, start=1, stop=3)
ppm <- 35
ms1_intensity_cutoff <- 100	          #approx. 0.01%

# General variables
mzml_files_pos <- NULL
mzml_names_pos <- NULL
mzml_times_pos <- NULL

# Load files
mzml_files_pos <- list.files(mzml_dir, pattern="*.mzML", recursive=T, full.names=T)
mzml_files_pos <- mzml_files_pos[grep(pol, mzml_files_pos, invert=FALSE)]
mzml_files_pos <- mzml_files_pos[grep("auto", mzml_files_pos, invert=FALSE)]

# Basenames of files without path and without extension
mzml_names_pos <- gsub('(.*)\\..*', '\\1', gsub('( |-|,)', '.', basename(mzml_files_pos)))

# Create phenodata based on species
mzml_pheno_pos <- data.frame(sample_name=mzml_names_pos, sample_group=gsub(pattern="\\..*", replacement="", x=mzml_names_pos, perl=TRUE))
mzml_pheno_samples_pos <- as.factor(mzml_pheno_pos$sample_group)
mzml_pheno_colors_pos <- data.frame(cols=c("steelblue2","steelblue3","orchid2","orchid3","red","palegreen2","palegreen3","brown2","brown3","blue","cyan3","cyan4"))
rownames(mzml_pheno_colors_pos) <- unique(mzml_pheno_pos$sample_group)
mzml_pheno_colors_samples_pos <- sapply(mzml_pheno_pos$sample_group, function(x) { x <- mzml_pheno_colors_pos[,1][which(x==unique(mzml_pheno_pos$sample_group))] } )

# Save timestamps of samples
for (i in 1:length(mzml_files_pos)) {
	fl <- mzR::openMSfile(mzml_files_pos[i])
	run_info <- mzR::runInfo(fl)
	mzR::close(fl)
	mzml_times_pos <- c(mzml_times_pos, run_info$startTimeStamp)
}

# Display MSn levels
mzml_msn_pos <- NULL
for (i in 1:length(mzml_files_pos)) {
	mzml_data_pos <- readMSData(mzml_files_pos[i], mode="onDisk")
	mzml_msn_pos <- rbind(mzml_msn_pos, t(as.matrix(table(msLevel(mzml_data_pos)))))
}
colnames(mzml_msn_pos) <- c("MS1", "MS2")
rownames(mzml_msn_pos) <- mzml_names_pos

# Plot MSn levels
pdf(file="plots/pos_msn_levels.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=16, family="Helvetica")
par(mfrow=c(2,1), mar=c(16,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
boxplot(mzml_msn_pos, main="Number of spectra")

model_boxplot <- boxplot(t(mzml_msn_pos[,2]), main="Number of MS2 spectra per sample", xaxt="n")
tick <- seq_along(model_boxplot$names)
axis(1, at=tick, labels=F)
text(tick, par("usr")[3]-par("usr")[3]/10, model_boxplot$names, adj=0, srt=270, xpd=T)
dev.off()



# ---------- Peak detection ----------
# Import raw data as MSnbase object
raw_data_pos <- readMSData(files=mzml_files_pos, pdata=new("NAnnotatedDataFrame", mzml_pheno_pos), mode="onDisk", centroided=TRUE)
table(msLevel(raw_data_pos))
head(fData(raw_data_pos)[, c("isolationWindowTargetMZ", "isolationWindowLowerOffset",
						     "isolationWindowUpperOffset", "msLevel", "retentionTime")])
write.csv(fData(raw_data_pos), file="data/pos_raw_data.csv", row.names=FALSE)

# Restrict data to 1020 seconds (17 minutes)
raw_data_pos <- filterRt(raw_data_pos, c(0, 1020))

# Inspect mz values per file
raw_mz_pos <- mz(raw_data_pos)
raw_mz_pos <- split(raw_mz_pos, f = fromFile(raw_data_pos))
print(length(raw_mz_pos))

# Get base peak chromatograms
chromas_pos <- chromatogram(raw_data_pos, aggregationFun="max")

# Plot chromatograms based on phenodata groups
pdf(file="plots/pos_chromas.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromas_pos, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col=mzml_pheno_colors_samples_pos)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=mzml_pheno_colors_pos[,1], legend=mzml_pheno_samples_pos)
dev.off()

# Get TICs
pdf(file="plots/pos_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(5,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
tics_pos <- split(tic(raw_data_pos), f=fromFile(raw_data_pos))
boxplot(tics_pos, names=mzml_names_pos, las=2, col=mzml_pheno_colors_samples_pos, ylab="intensity", main="Total ion current")
dev.off()

# Grouping/binning the samples based on similarity of their base peak chromatogram to spot potentially problematic samples
chromas_bin_pos <- MSnbase::bin(chromas_pos, binSize=2)
chromas_bin_cor_pos <- cor(log2(do.call(cbind, lapply(chromas_bin_pos, intensity))))
colnames(chromas_bin_cor_pos) <- rownames(chromas_bin_cor_pos) <- raw_data_pos$sample_name
pdf(file="plots/pos_chromas_bin_cor.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_cor_pos)
dev.off()

# Assess retention times and intensities of first file
head(rtime(chromas_pos[1, 1]))
head(intensity(chromas_pos[1, 1]))

# Inspect peak width of standard compound for defining base peakwidth parameter below
pdf(file="plots/pos_chromas_standard_peakwidth_348-350.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromatogram(raw_data_pos, mz=c(215, 217), rt=c(214, 219)), col=mzml_pheno_colors_samples_pos, main = "EIC of Kinetin")
plot(chromatogram(raw_data_pos, mz=c(284, 286), rt=c(546, 553)), col=mzml_pheno_colors_samples_pos, main = "EIC of Biochanin A")
plot(chromatogram(raw_data_pos, mz=c(246, 248), rt=c(338, 346)), col=mzml_pheno_colors_samples_pos, main = "EIC of N-(3-Indolylacetyl)-L-Ala")

plot(chromatogram(raw_data_pos, mz=c(280, 282), rt=c(670, 676)), col=mzml_pheno_colors_samples_pos, main = "EIC of Radulanin A")
plot(chromatogram(raw_data_pos, mz=c(324, 326), rt=c(676, 682)), col=mzml_pheno_colors_samples_pos, main = "EIC of Radulanin H")
plot(chromatogram(raw_data_pos, mz=c(296, 298), rt=c(589, 597)), col=mzml_pheno_colors_samples_pos, main = "EIC of Radulanin L")

plot(chromatogram(raw_data_pos, mz=c(282, 284), rt=c(615, 621)), col=mzml_pheno_colors_samples_pos, main = "EIC of 4-prenyldihydropinosylvin")
plot(chromatogram(raw_data_pos, mz=c(350, 352), rt=c(754, 761)), col=mzml_pheno_colors_samples_pos, main = "EIC of Bibenzyl-CBG")
plot(chromatogram(raw_data_pos, mz=c(340, 342), rt=c(745, 751)), col=mzml_pheno_colors_samples_pos, main = "EIC of Amorfrutin A")

plot(chromatogram(raw_data_pos, mz=c(348, 350), rt=c(0, 1200)), col=mzml_pheno_colors_samples_pos, main = "EIC of 348-350")
dev.off()

# Peak detection in MS1 data
if (polarity=="negative") {
	ms1_params_neg <- CentWaveParam(ppm=13, mzCenterFun="mean", peakwidth=c(4, 33), prefilter=c(4, 200), mzdiff=0.0023, snthresh=11, noise=0, integrate=1,
								    firstBaselineCheck=TRUE, verboseColumns=TRUE, fitgauss=FALSE, roiList=list(), roiScales=numeric())
} else {
	ms1_params_pos <- CentWaveParam(ppm=32, mzCenterFun="mean", peakwidth=c(4, 32), prefilter=c(2, 100), mzdiff=0.0111, snthresh=8, noise=0, integrate=1,
							        firstBaselineCheck=TRUE, verboseColumns=FALSE, fitgauss=FALSE, roiList=list(), roiScales=numeric())
}
ms1_data_pos <- findChromPeaks(raw_data_pos, param=ms1_params_pos)

# Per file summary
ms1_summary_pos <- lapply(split.data.frame(chromPeaks(ms1_data_pos), f=chromPeaks(ms1_data_pos)[, "sample"]), FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )
ms1_summary_pos <- do.call(rbind, ms1_summary_pos)
rownames(ms1_summary_pos) <- basename(fileNames(ms1_data_pos))
print(ms1_summary_pos)
table(msLevel(ms1_data_pos))
write.csv(as.data.frame(table(msLevel(ms1_data_pos))), file="data/pos_ms1_data.csv", row.names=FALSE)

# To get a global overview of the peak detection we can plot the frequency of identified peaks per file along the retention time axis. This allows to identify time periods along the MS run in which a higher number of peaks was identified and evaluate whether this is consistent across files.
pdf(file="plots/pos_ms1_data.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,18,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plotChromPeakImage(ms1_data_pos, main="Frequency of identified peaks per RT")
dev.off()

# Group peaks
if (polarity=="negative") {
	ms1_data_neg <- groupChromPeaks(ms1_data_neg, param=PeakDensityParam(sampleGroups=ms1_data_neg$sample_group, minFraction=0.7, bw=2.5))
} else {
	ms1_data_pos <- groupChromPeaks(ms1_data_pos, param=PeakDensityParam(sampleGroups=ms1_data_pos$sample_group, minFraction=0.7, bw=22))
}

# RT correction
if (polarity=="negative") {
	ms1_data_neg <- adjustRtime(ms1_data_neg, param=PeakGroupsParam(minFraction=0.7,smooth="loess",span=0.2,family="gaussian"))
} else {
	ms1_data_pos <- adjustRtime(ms1_data_pos, param=PeakGroupsParam(minFraction=0.7,smooth="loess",span=0.2,family="gaussian"))
}

# Plot the difference of raw and adjusted retention times
pdf(file="plots/pos_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(chromas_pos, peakType="none", main="Raw chromatograms")
plotAdjustedRtime(ms1_data_pos, lwd=2, main="Retention Time correction")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
dev.off()

# Group peaks
if (polarity=="negative") {
	ms1_data_neg <- groupChromPeaks(ms1_data_neg, param=PeakDensityParam(sampleGroups=ms1_data_neg$sample_group, minFraction=0.7, bw=2.5))
} else {
	ms1_data_pos <- groupChromPeaks(ms1_data_pos, param=PeakDensityParam(sampleGroups=ms1_data_pos$sample_group, minFraction=0.7, bw=22))
}

# Get integrated peak intensity per feature/sample
print(head(featureValues(ms1_data_pos, value="into")))

# Fill peaks
ms1_data_pos <- fillChromPeaks(ms1_data_pos, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
head(featureValues(ms1_data_pos))
head(featureSummary(ms1_data_pos, group=ms1_data_pos$sample_group))

# Evaluate grouping
pdf(file="plots/pos_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca_pos <- prcomp(t(na.omit(log2(featureValues(ms1_data_pos, value="into")))), center=TRUE)
plot(ms1_pca_pos$x[, 1], ms1_pca_pos$x[,2], pch=19, main="PCA: Grouping of samples",
	 xlab=paste0("PC1: ", format(summary(ms1_pca_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(ms1_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 1], ms1_pca_pos$x[,2], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()

# Show peaks
tail(chromPeaks(ms1_data_pos))
tail(chromPeakData(ms1_data_pos))

# Show process history
processHistory(ms1_data_pos)



# ---------- Build MS1 feature tables ----------
# Build feature matrix
ms1_matrix_pos <- featureValues(ms1_data_pos, method="medret", value="into")
colnames(ms1_matrix_pos) <- mzml_names_pos
dim(ms1_matrix_pos)
feat_list_pos <- t(ms1_matrix_pos)

# Build feature summary
ms1_summary_pos <- featureSummary(ms1_data_pos)
ms1_def_pos <- featureDefinitions(ms1_data_pos)

# Missing value imputation
feat_list_pos[is.na(feat_list_pos)] <- median(na.omit(as.numeric(unlist(feat_list_pos))))

# Transform data
feat_list_pos <- log2(feat_list_pos)

# Missing value imputation
feat_list_pos[which(is.na(feat_list_pos))] <- median(na.omit(as.numeric(unlist(feat_list_pos))))

#feat_list_pos <- scale(feat_list_pos, scale=FALSE, center=FALSE)
#feat_list_pos[is.na(feat_list_pos)] <- 0
#feat_list_pos[which(feat_list_pos < 0)] <- 0
#feat_list_pos[is.infinite(feat_list_pos)] <- 0
#feat_list_pos <- feat_list_pos[!apply(feat_list_pos, MARGIN=1, function(x) max(x,na.rm=TRUE) == min(x,na.rm=TRUE)),]

# Plot histogram
pdf(file="plots/pos_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(feat_list_pos), main="Histogram of feature table")
dev.off()

# PCA
pdf(file="plots/pos_ms1_feature_table_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca_pos <- prcomp(feat_list_pos, center=TRUE)
plot(ms1_pca_pos$x[, 1], ms1_pca_pos$x[,2], pch=19, main="PCA of feature table",
	 xlab=paste0("PC1: ", format(summary(ms1_pca_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(ms1_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 1], ms1_pca_pos$x[,2], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()

# Create single 0/1 matrix
bina_list_pos <- t(ms1_matrix_pos)
bina_list_pos[is.na(bina_list_pos)] <- 1
bina_list_pos <- log2(bina_list_pos)
bina_list_pos[bina_list_pos < log2(ms1_intensity_cutoff)] <- 0
bina_list_pos[bina_list_pos != 0] <- 1

# Only unique compounds in group mzml_pheno$ and not the others
uniq_list_pos <- apply(X=bina_list_pos, MARGIN=2, FUN=function(x) { if (length(unique(mzml_pheno_pos$sample_group[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
colnames(uniq_list_pos) <- colnames(bina_list_pos)
rownames(uniq_list_pos) <- rownames(bina_list_pos)

# Create data frame
model_div_pos             <- data.frame(features=apply(X=bina_list_pos, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div_pos$richness    <- apply(X=bina_list_pos, MARGIN=1, FUN=function(x) { sum(x) } )
#model_div_pos$menhinick   <- apply(X=bina_list_pos, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div_pos$shannon     <- apply(X=feat_list_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div_pos$pielou      <- apply(X=scale(feat_list_pos, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div_pos$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list_pos, species), index="chao")
model_div_pos$simpson     <- apply(X=feat_list_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div_pos$inverse     <- apply(X=feat_list_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div_pos$fisher      <- apply(X=feat_list_pos, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div_pos$unique      <- apply(X=uniq_list_pos, MARGIN=1, FUN=function(x) { sum(x) })

# Remove NAs if present
model_div_pos[is.na(model_div_pos)] <- 0



# ############################## MS2 ##############################



# ---------- MS2 spectra detection ----------
# Estimate precursor intensity
precursor_intensity_pos <- estimatePrecursorIntensity(ms1_data_pos)
print(head(na.omit(precursor_intensity_pos)))

# Reconstruct MS2 spectra from MS1 data
ms2_data_pos <- chromPeakSpectra(ms1_data_pos, msLevel=2L, return.type="Spectra")
print(ms2_data_pos)
print(length(ms2_data_pos$peak_id))

# Extract all usable MS2 spectra
ms2_spectra_pos <- list()
#for (i in 1:nrow(ms1_def_pos)) {
ms2_spectra_pos <- foreach(i=1:nrow(ms1_def_pos)) %dopar% {
	# Extract existing MS2 spectra for feature
	feature_of_interest <- ms1_def_pos[i, "mzmed"]
	peaks_of_interest <- chromPeaks(ms1_data_pos, mz=feature_of_interest, ppm=ppm)
	if (length(which(ms2_data_pos$peak_id %in% rownames(peaks_of_interest))) > 0) {
		spectra_of_interest <- ms2_data_pos[ms2_data_pos$peak_id %in% rownames(peaks_of_interest)]
		combined_spectra_of_interest <- filterIntensity(spectra_of_interest, intensity=c(1,Inf), backend=MsBackendDataFrame)
		combined_spectra_of_interest <- setBackend(combined_spectra_of_interest, backend=MsBackendDataFrame())
		combined_spectra_of_interest <- Spectra::combineSpectra(combined_spectra_of_interest, FUN=combinePeaks, ppm=ppm, peaks="union", minProp=0.8, intensityFun=median, mzFun=median, backend=MsBackendDataFrame)#f=rownames(peaks_of_interest))
		combined_spectra_peaks <- as.data.frame(peaksData(combined_spectra_of_interest)[[1]])
		#Spectra::plotSpectra(combined_spectra_of_interest)
		#plot(x=combined_spectra_peaks[,1], y=combined_spectra_peaks[,2], type="h", xlab="m/z", ylab="intensity", main=paste("Precursor m/z",combined_spectra_of_interest@backend@spectraData$precursorMz[[1]]))
		#length(spectra_of_interest$peak_id)
		
		#ms2_spectra_pos[[i]] <- combined_spectra_of_interest
		return(combined_spectra_of_interest)
	}
}

# Remove empty spectra
names(ms2_spectra_pos) <- rownames(ms1_def_pos)
ms2_spectra_pos <- ms2_spectra_pos[lengths(ms2_spectra_pos) != 0]

# Relate all MS2 spectra to MS1 precursors
ms1_def_pos$has_ms2 <- as.integer(rownames(ms1_def_pos) %in% names(ms2_spectra_pos))



# ---------- Annotate MS2 spectra ----------
# Save all MS2 spectra in MSP file
msp_text <- NULL
alignment_id <- 0
for (i in names(ms2_spectra_pos)) {
	alignment_id <- alignment_id + 1
	
	msp_text <- c(msp_text, paste("NAME:", i))
	msp_text <- c(msp_text, paste("AlignmentID:", alignment_id))
	msp_text <- c(msp_text, paste("RETENTIONTIME:", ms1_def_pos[i, "rtmed"]))
	msp_text <- c(msp_text, paste("PRECURSORMZ:", ms1_def_pos[i, "mzmed"]))
	msp_text <- c(msp_text, paste("METABOLITENAME:", i))
	if (polarity == "positive") {
		msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
	} else {
		msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
	}
	msp_text <- c(msp_text, paste("NumPeaks:", nrow(as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]]))))
	msp_text <- c(msp_text, paste(as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$mz, as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$intensity, sep="\t"))
	msp_text <- c(msp_text, "")
	
}

# Write MSP file
cat(msp_text, file="data/pos_ms2_spectra.msp", sep="\n")


# Save all MS2 spectra in MGF file
mgf_text <- NULL
for (i in names(ms2_spectra_pos)) {
	mgf_text <- c(mgf_text, paste0("COM=", i))
	mgf_text <- c(mgf_text, "BEGIN IONS")
	mgf_text <- c(mgf_text, "MSLEVEL=2")
	mgf_text <- c(mgf_text, paste0("TITLE=", i))
	mgf_text <- c(mgf_text, paste0("RTINSECONDS=", ms1_def_pos[i, "rtmed"]))
	mgf_text <- c(mgf_text, paste0("PEPMASS=", ms1_def_pos[i, "mzmed"]))
	if (polarity == "positive") {
		mgf_text <- c(mgf_text, paste0("CHARGE=", "1+"))
	} else {
		mgf_text <- c(mgf_text, paste0("CHARGE=", "1-"))
	}
	mgf_text <- c(mgf_text, paste(as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$mz, as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$intensity, sep=" "))
	mgf_text <- c(mgf_text, "END IONS")
	mgf_text <- c(mgf_text, "")
}

# Write MGF file
cat(mgf_text, file="data/pos_ms2_spectra.mgf", sep="\n")



# ---------- Classify MS2 spectra with MetFamily ----------
# Load MetFamily
source("MetFamily/MetFamily.R")
source("MetFamily/R_packages.R")
source("MetFamily/FragmentMatrixFunctions.R")
source("MetFamily/DataProcessing.R")
source("MetFamily/Annotation.R")
source("MetFamily/Classifiers.R")

# Define MetFamily classifier
if (polarity == "positive") {
	classifier_name_pos <- "2019-04-17_MSMS_pos_21908_MoNA_Classifier"
} else {
	classifier_name_neg <- "2019-06-14_neg_35009_MoNA_Classifier"
}

# Read CHEMONT ontology
obo <- ontologyIndex::get_ontology(file=paste0(mzml_dir, "ChemOnt_2_1.obo"), extract_tags="minimal")

# Associate CHEMONTIDs with classes names
classes_names <- obo$name
classes_id <- NULL
for (i in 1:length(obo$id)) {
	classes_id <- c(classes_id, as.character(obo$id[which(as.character(obo$name) == gsub(x=obo$name[i], pattern=".*; ", replacement=""))]))
}

# Apply MetFamily classifier to all spectra
classifier_metfamily_pos <- tryCatch(
	obj <- applyClassifierMs2(classifierFile = paste0("MetFamily/", classifier_name_pos, ".RData"),
							  propertiesFile = paste0("MetFamily/", classifier_name_pos, ".txt"),
							  fileMs1Path = NULL,
							  fileMs2Path = "data/pos_ms2_spectra.msp",
							  fileClasses = "raw/massbank_classes.txt",
							  minimumIntensityOfMaximalMS2peak = msms.intensity.threshold,
							  minimumProportionOfMS2peaks = 0.05,
							  mzDeviationAbsolute_grouping = mzabs,
							  mzDeviationInPPM_grouping = mzppm),
	error = function(c) NULL,
	warning = function(c) obj,
	message = function(c) obj)

# Save classified objects
for (i in 1:1) {
	write.table(x=classifier_metfamily_pos, file=paste0("data/pos_classifier_metfamily.tsv"), sep="\t", quote=TRUE, row.names=FALSE, dec=".")
}

# Re-Load classifier
load(paste0("MetFamily/", classifier_name_pos, ".RData"))

# Read classifiers
classifier_metfamily_pos <- list()
for (i in 1:1) {
	classifier_metfamily_pos[[i]] <- read.table(file=paste0("data/pos_classifier_metfamily.tsv"), header=TRUE, sep="\t", quote="\"", fill=FALSE, dec=".", stringsAsFactors=FALSE)
	
	# Sort data frame with increasing p-value
	classifier_metfamily_pos[[i]] <- classifier_metfamily_pos[[i]][order(classifier_metfamily_pos[[i]]$pValue,decreasing=FALSE), ]
}

# Apply classified classes onto feature table
ms1_def_pos$primary_class <- ""
ms1_def_pos$alternative_classes <- ""

for (i in 1:1) {
	for (j in unique(classifier_metfamily_pos[[i]]$Metabolite.name)) {
		obj <- classifier_metfamily_pos[[i]][classifier_metfamily_pos[[i]]$Metabolite.name %in% j, "Annotation..putative."]
		ms1_def_pos[j, "primary_class"] <- obj[1]
		ms1_def_pos[j, "alternative_classes"] <- paste0(obj[2:length(obj)], collapse="|")
	}
}

# Read SIRIUS/CANOPUS classifier
if (CANOPUS == TRUE) {
	classifier_canopus_pos <- read.table(file=paste0("data/pos_ms2_canopus_summary.tsv"), header=TRUE, sep="\t", quote="\"", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
	classifier_canopus_pos$Metabolite.name <- gsub(x=classifier_canopus_pos$name, pattern=".*_", replacement="")
	classifier_canopus_pos$primary_class <- paste("Organic compounds", classifier_canopus_pos$superclass, classifier_canopus_pos$class, classifier_canopus_pos$subclass, classifier_canopus_pos$level.5, sep="; ")
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_pos$Annotation..putative. <- classifier_canopus_pos$primary_class
	classifier_canopus_pos$alternative_classes <- classifier_canopus_pos$all.classifications
	
	for (j in unique(classifier_canopus_pos$Metabolite.name)) {
		obj <- classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% j, "Annotation..putative."]
		ms1_def_pos[j, "primary_class"] <- obj[1]
	}
}



# ---------- Annotate MS2 spectra for samples ----------
# Write MSP for all samples
foreach (j=mzml_names_pos) %dopar% {
#for (j in mzml_names_pos) {
	# Iterate through all MS2 spectra found in the sample
	msp_text <- NULL
	alignment_id <- 0
	for (i in names(ms2_spectra_pos)) {
		if (bina_list_pos[j,i] != 0) {
			alignment_id <- alignment_id + 1
			
			msp_text <- c(msp_text, paste("NAME:", i))
			msp_text <- c(msp_text, paste("AlignmentID:", alignment_id))
			msp_text <- c(msp_text, paste("RETENTIONTIME:", ms1_def_pos[i, "rtmed"]))
			msp_text <- c(msp_text, paste("PRECURSORMZ:", ms1_def_pos[i, "mzmed"]))
			msp_text <- c(msp_text, paste("METABOLITENAME:", i))
			msp_text <- c(msp_text, paste("COMPOUND_CLASS:", ms1_def_pos[i, "primary_class"]))
			msp_text <- c(msp_text, paste("ALTERNATIVE_CLASSES:", ms1_def_pos[i, "alternative_classes"]))
			if (polarity == "positive") {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
			} else {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
			}
			msp_text <- c(msp_text, paste("NumPeaks:", nrow(as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]]))))
			msp_text <- c(msp_text, paste(as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$mz, as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$intensity, sep="\t"))
			msp_text <- c(msp_text, "")
		}
	}
	
	# Write MSP file
	cat(msp_text, file=paste0("msp/",j,".msp"), sep="\n")
}



# ---------- Diversity of MS2 classes ----------
# Create classifier object for each sample
#classifiers_pos <- list()
#for (i in mzml_names_pos) {
#	obj <- names(which(bina_list_pos[rownames(bina_list_pos)==i, colnames(bina_list_pos) %in% classifier_metfamily_pos[[1]]$Metabolite.name] > 0))
#	classifiers_pos[[i]] <- classifier_metfamily_pos[[1]][classifier_metfamily_pos[[1]]$Metabolite.name %in% obj, ]
#}

# Create CANOPUS classifier object for each sample
if (CANOPUS == TRUE) {
	classifiers_pos <- list()
	for (i in mzml_names_pos) {
		obj <- names(which(bina_list_pos[rownames(bina_list_pos)==i, colnames(bina_list_pos) %in% classifier_canopus_pos$Metabolite.name] > 0))
		classifiers_pos[[i]] <- classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% obj, ]
	}
}

# Diversity of classes per sample
div_classes_samples_pos <- NULL
for (i in mzml_names_pos) {
	obj <- table(classifiers_pos[[i]][,"Annotation..putative."])
	obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
	if (is.null(div_classes_samples_pos)) {
		div_classes_samples_pos <- obj
	} else {
		div_classes_samples_pos <- merge(div_classes_samples_pos, obj, by="classes", all.x=TRUE, all.y=TRUE)
	}
}
rownames(div_classes_samples_pos) <- div_classes_samples_pos$classes
div_classes_samples_pos <- div_classes_samples_pos[, -which(colnames(div_classes_samples_pos)=="classes")]
colnames(div_classes_samples_pos) <- mzml_names_pos

# Diversity of classes
div_classes_pos <- div_classes_samples_pos
div_classes_pos[is.na(div_classes_pos)] <- 0
div_classes_pos <- apply(X=div_classes_pos, MARGIN=1, FUN=function(x) { sum(x) })
div_classes_pos <- data.frame(row.names=names(div_classes_pos), frequency=as.numeric(div_classes_pos))

# Plot diversity of classes in all samples
pdf(file="plots/pos_ms2_classes_diversity.pdf", encoding="ISOLatin1", pointsize=8, width=6, height=14, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,15,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(div_classes_pos$frequency, names.arg=gsub('.*; ','',rownames(div_classes_pos)), las=1, horiz=TRUE, xlab="frequency", main="Diversity of compound classes", col=rainbow(n=nrow(div_classes_pos), alpha=0.6))
dev.off()

# Sunburst plot of classes of all samples
pdf(file="plots/pos_ms2_classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(div_classes_pos), div_classes_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_classes_pos), "Number of spectra"=div_classes_pos$frequency), file="plots/pos_ms2_classes_sunburst.csv", row.names=FALSE)

# Classes
classes_order_pos <- readLines(con="raw/massbank_classes.txt")
classes_order_pos <- classes_order_pos[which(grepl(pattern="^Organic compounds", x=classes_order_pos))]
classes_order_pos <- gsub(x=classes_order_pos, pattern=":", replacement="")
classes_order_pos <- gsub(x=classes_order_pos, pattern="/", replacement="; ")

classes_pos <- rownames(div_classes_pos)
classes_pos <- classes_pos[which(grepl(pattern="^Organic compounds", x=classes_pos))]
classes_pos <- gsub(x=classes_pos, pattern=":", replacement="")
classes_pos <- gsub(x=classes_pos, pattern="/", replacement="; ")

classes_pos <- classes_order_pos[which(classes_order_pos %in% classes_pos)]

# Determine how many spectra were classified
spectra_number_pos <- ncol(feat_list_pos)
spectra_classified_pos <- sum(div_classes_pos)/length(mzml_names_pos)
print(paste("Number of spectra:", spectra_number_pos))
print(paste("Number of spectra classified:", spectra_classified_pos))
print(paste("Number of unclassified spectra:", spectra_number_pos - spectra_classified_pos))
print(paste("Number of classes:", length(classes_order_pos)))
print(paste("Number of classes with entities:", length(classes_pos)))
print(paste("Number of classes without entities:", length(classes_order_pos) - length(classes_pos)))

# Classifiers
classifiers_class_pos <- get(load(paste0("MetFamily/", classifier_name_pos, ".RData")))

# Area Under Curve
classifier_auc_pos <- as.numeric(unlist(lapply(X=classes_pos, FUN = function(x) { x <- classifiers_class_pos[[x]]$AUC } )))

# Area under Precision Recall Curve
classifier_aucpr_pos <- as.numeric(unlist(lapply(X=classes_pos, FUN = function(x) { x <- classifiers_class_pos[[x]]$AUC_PR } )))

# True Positive Rate for False Negative Rate of 5 Percent
classifier_fnr_pos <- as.numeric(unlist(lapply(X=classes_pos, FUN = function(x) { x <- classifiers_class_pos[[x]]$TNR_for_FNR_of_5Percent } )))

# True Positive Rate for False Positive Rate of 5 Percent
classifier_fpr_pos <- as.numeric(unlist(lapply(X=classes_pos, FUN = function(x) { x <- classifiers_class_pos[[x]]$TPR_for_FPR_of_5Percent } )))

# Save table with AUC-PR and TPR-FPR rates
write.table(x=data.frame(compound_class=classes_pos, "AUC"=classifier_auc_pos, "AUC-PR"=classifier_aucpr_pos, "TPR-FPR"=classifier_fnr_pos, "TPR-FPR"=classifier_fpr_pos),
			file="data/pos_ms2_classifier_performance.csv",
			sep=";", quote=TRUE, row.names=FALSE, dec=".")

# Remove object to save memory
rm(classifiers_class_pos)



# ---------- Build diversity objects ----------
# Imputation of NA with zeros
div_classes_pos[is.na(div_classes_pos)] <- 0
div_classes_samples_pos[is.na(div_classes_samples_pos)] <- 0

# Classification list for statistics
class_list_pos <- as.data.frame(t(div_classes_samples_pos))
class_list_pos[is.na(class_list_pos)] <- 0

# Log Transform
#class_list_pos <- log2(class_list_pos + 1)

# Only keep class names
colnames(class_list_pos) <- gsub(x=colnames(class_list_pos), pattern='.*; ', replacement='')



# ---------- Classification at level of superclasses ----------
# Use CANOPUS
if (CANOPUS == TRUE) {
	#classifier_canopus_pos$class,
	classifier_canopus_pos$primary_class <- paste("Organic compounds", classifier_canopus_pos$superclass, classifier_canopus_pos$most.specific.class, sep="; ")
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_pos$Annotation..putative. <- classifier_canopus_pos$primary_class
	
	# Create CANOPUS classifier object for each sample
	superclassifiers_pos <- list()
	for (i in mzml_names_pos) {
		obj <- names(which(bina_list_pos[rownames(bina_list_pos)==i, colnames(bina_list_pos) %in% classifier_canopus_pos$Metabolite.name] > 0))
		superclassifiers_pos[[i]] <- classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% obj, ]
	}
	
	# Diversity of classes per sample
	div_superclasses_samples_pos <- NULL
	for (i in mzml_names_pos) {
		obj <- table(superclassifiers_pos[[i]][,"Annotation..putative."])
		obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
		if (is.null(div_superclasses_samples_pos)) {
			div_superclasses_samples_pos <- obj
		} else {
			div_superclasses_samples_pos <- merge(div_superclasses_samples_pos, obj, by="classes", all.x=TRUE, all.y=TRUE)
		}
	}
	rownames(div_superclasses_samples_pos) <- div_superclasses_samples_pos$classes
	div_superclasses_samples_pos <- div_superclasses_samples_pos[, -which(colnames(div_superclasses_samples_pos)=="classes")]
	colnames(div_superclasses_samples_pos) <- mzml_names_pos
} else {
	# Make superclasses at level 2
	superclass_level_pos <- 2
	div_superclasses_samples_names_pos <- NULL
	for (i in c(1:superclass_level_pos)) {
		div_superclasses_samples_names_pos <- c(div_superclasses_samples_names_pos, lapply(X=strsplit(rownames(div_classes_samples_pos), '; '), FUN=function(x) { gsub(x=paste(x[1:i],sep='',collapse='; '),pattern='; NA',replacement='') }))
	}
	div_superclasses_samples_pos <- data.frame()
	for (i in c(1:ncol(div_classes_samples_pos))) div_superclasses_samples_pos <- rbind(div_superclasses_samples_pos, rep(0, length(unique(div_superclasses_samples_names_pos))))
	div_superclasses_samples_pos <- t(div_superclasses_samples_pos)
	colnames(div_superclasses_samples_pos) <- colnames(div_classes_samples_pos)
	rownames(div_superclasses_samples_pos) <- unique(div_superclasses_samples_names_pos)
	for (i in rownames(div_superclasses_samples_pos)) {
		for (j in c(1:ncol(div_classes_samples_pos))) {
			div_superclasses_samples_pos[rownames(div_superclasses_samples_pos)==i, j] <- sum(div_classes_samples_pos[grep(x=rownames(div_classes_samples_pos), pattern=i), j])
		}
	}
}

# Diversity of superclasses
div_superclasses_pos <- div_superclasses_samples_pos
div_superclasses_pos[is.na(div_superclasses_pos)] <- 0
div_superclasses_pos <- apply(X=div_superclasses_pos, MARGIN=1, FUN=function(x) { sum(x) })
div_superclasses_pos <- data.frame(row.names=names(div_superclasses_pos), frequency=as.numeric(div_superclasses_pos))

# Sunburst plot
pdf(file="plots/pos_ms2_superclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_superclasses_pos), div_superclasses_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_superclasses_pos), "Number of spectra"=div_superclasses_pos$frequency), file="plots/pos_ms2_superclasses_sunburst.csv", row.names=FALSE)

# Imputation of NA with zeros
div_superclasses_pos[is.na(div_superclasses_pos)] <- 0
div_superclasses_samples_pos[is.na(div_superclasses_samples_pos)] <- 0

# Classification list for statistics
superclass_list_pos <- as.data.frame(t(div_superclasses_samples_pos))
superclass_list_pos[is.na(superclass_list_pos)] <- 0

# Log Transform
#superclass_list_pos <- log2(superclass_list_pos + 1)

# Only keep superclass names
colnames(superclass_list_pos) <- gsub(x=colnames(superclass_list_pos), pattern='.*; ', replacement='')



# ############################## MS1 statistics ##############################



# ---------- Histogram ----------
pdf(file="plots/pos_ms1_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(feat_list_pos))
dev.off()



# ---------- Variation partitioning ----------
mzml_pheno_growth_stress_samples_pos <- as.factor(c("Stress", "Stress", "Stress", "Stress", "Stress", "Stress", 
													"Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth", "Growth",
													"Stress", "Stress", "Stress", "Stress", "Stress", "Stress",
													"Growth", "Growth", "Growth", "Growth", "Growth", "Growth",
													"Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress", "Stress"))
model_varpart_pos <- varpart(scale(feat_list_pos), ~ mzml_pheno_samples_pos, ~ mzml_pheno_growth_stress_samples_pos)

# Plot results
pdf(file="plots/pos_ms1_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart_pos, Xnames=c("samples","hormone treatment"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- Diversity ----------
# Plot unique features
pdf(paste("plots/pos_ms1_richness_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div_pos$unique ~ mzml_pheno_samples_pos, col=mzml_pheno_colors_pos$cols, names=NA, main="Number of unique features", xlab="treatment", ylab="number of unique features")
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_pos$sample_group), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div_pos$unique, term=as.factor(mzml_pheno_samples_pos))
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot total features
pdf(paste("plots/pos_ms1_richness_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div_pos$features ~ mzml_pheno_samples_pos, col=mzml_pheno_colors_pos$cols, names=NA, main="Number of features", xlab="treatment", ylab="number of features")
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_pos$sample_group), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div_pos$features, term=as.factor(mzml_pheno_samples_pos))
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Shannon index
pdf(paste("plots/pos_ms1_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div_pos$shannon ~ mzml_pheno_samples_pos, col=mzml_pheno_colors_pos$cols, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_pos$sample_group), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div_pos$shannon, term=as.factor(mzml_pheno_samples_pos))
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Pielou evenness for species
pdf(paste("plots/pos_ms1_diversity_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div_pos$pielou ~ mzml_pheno_samples_pos, col=mzml_pheno_colors_pos$cols, names=NA, main="Pielou\'s evenness", xlab="treatment", ylab="Pielou diversity index (J)")
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=unique(mzml_pheno_pos$sample_group), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div_pos$pielou, term=as.factor(mzml_pheno_samples_pos))
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()



# ---------- Variable Selection ----------
# Random Forest
sel_rf_pos <- f.select_features_random_forest(feat_matrix=feat_list_pos, sel_factor=as.factor(mzml_pheno_samples_pos), sel_colors=mzml_pheno_colors_pos$cols, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/pos_ms1_select_rf_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_rf_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=feat_list_pos, sel_feat=sel_rf_pos$`_selected_variables_`, filename="plots/pos_ms1_select_rf.pdf", main="Random Forest")

# PLS
sel_pls_pos <- f.select_features_pls(feat_matrix=feat_list_pos, sel_factor=as.factor(mzml_pheno_samples_pos), sel_colors=mzml_pheno_colors_pos$cols, components=2, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/pos_ms1_select_pls_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_pls_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=feat_list_pos, sel_feat=sel_pls_pos$`_selected_variables_`, filename="plots/pos_ms1_select_pls.pdf", main="PLS")

# sPLS-DA
# First: Variation partitioning
model_varpart_pos <- varpart(scale(feat_list_pos), ~ mzml_pheno_samples_pos, ~ mzml_pheno_samples_pos)

# 10% of features correlate with factor
model_varpart_corr_pos <- trunc(model_varpart_pos$part$indfract$Adj.R.squared[2] * ncol(feat_list_pos))
model_splsda_keepx_pos <- trunc(seq(model_varpart_corr_pos / length(unique(mzml_pheno_samples_pos)), model_varpart_corr_pos / length(unique(mzml_pheno_samples_pos))^2,length.out=length(unique(mzml_pheno_samples_pos))))

sel_splsda_pos <- f.select_features_splsda(feat_matrix=feat_list_pos, sel_colors=mzml_pheno_colors_pos$cols, sel_factor=mzml_pheno_samples_pos, tune_components=length(unique(mzml_pheno_samples_pos)) - 1, sel_components=c(3), folds_number=10, keepx=model_splsda_keepx_pos, plot_roc_filename="plots/pos_ms1_select_splsda_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_splsda_pos$'_selected_variables_')))
f.heatmap.selected_features(feat_list=feat_list_pos, sel_feat=sel_splsda_pos$'_selected_variables_', filename="plots/pos_ms1_select_splsda.pdf", main="sPLS-DA")



# ############################## MS2 Classes statistics ##############################



# ---------- Histogram ----------
pdf(file="plots/pos_ms2_classes_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(unlist(class_list_pos)), main="Histogram", xlab="class_list")
dev.off()



# ---------- Variation partitioning ----------
model_varpart_pos <- varpart(class_list_pos, ~ mzml_pheno_samples_pos, ~ mzml_pheno_samples_pos)

# Plot results
pdf(file="plots/pos_ms2_classes_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart_pos, Xnames=c("sample group","samples"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca_pos <- prcomp(class_list_pos)

pdf(paste("plots/pos_ms2_classes_pca.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca_pos$x[, 1], model_pca_pos$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca_pos)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, bg=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(model_pca_pos$x[, 1], model_pca_pos$x[,2], labels=mzml_names_pos, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()


# ---------- Feature selection of classes ----------
# Tukey on ANOVA
#sel_aov_pos <- f.select_features_aov(feat_class=mzml_pheno_samples_pos, feat_list=class_list_pos, conf_level=0.95, pvalue=0.99)
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_aov_pos)))
#f.heatmap.selected_features(feat_list=class_list_pos, sel_feat=sel_aov_pos, filename="plots/pos_ms2_classes_select_aov.pdf", main="Tukey on ANOVA")

# Linear Model
#sel_lm_pos <- f.select_features_linear_model(feat_matrix=class_list_pos, sel_factor=mzml_pheno_samples_pos, sel_formula=feature ~ mzml_pheno_samples_pos, pvalue=0.05)
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_lm_pos)))
#f.heatmap.selected_features(feat_list=class_list_pos, sel_feat=sel_lm_pos, filename="plots/pos_ms2_classes_select_lm.pdf", main="Linear Model")

# Random Forest
sel_rf_pos <- f.select_features_random_forest(feat_matrix=class_list_pos, sel_factor=mzml_pheno_samples_pos, sel_colors=mzml_pheno_colors_pos$cols, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/pos_ms2_classes_select_rf_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_rf_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=class_list_pos, sel_feat=sel_rf_pos$`_selected_variables_`, filename="plots/pos_ms2_classes_select_rf.pdf", main="Random Forest")

# PLS
sel_pls_pos <- f.select_features_pls(feat_matrix=class_list_pos, sel_factor=mzml_pheno_samples_pos, sel_colors=mzml_pheno_colors_pos$cols, components=2, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/pos_ms2_classes_select_pls_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_pls_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=class_list_pos, sel_feat=sel_pls_pos$`_selected_variables_`, filename="plots/pos_ms2_classes_select_pls.pdf", main="PLS")

# sPLS-DA
# First: Variation partitioning
#model_varpart_pos <- varpart(scale(class_list_pos), ~ mzml_pheno_samples_pos, ~ mzml_pheno_samples_pos)

# 10% of features correlate with factor
#model_varpart_corr_pos <- trunc(model_varpart_pos$part$indfract$Adj.R.squared[2] * ncol(class_list_pos))
#model_splsda_keepx_pos <- trunc(seq(model_varpart_corr_pos, 1,length.out=length(unique(mzml_pheno_samples_pos))))

#sel_splsda_pos <- f.select_features_splsda(feat_matrix=class_list_pos, sel_colors=mzml_pheno_colors_pos$cols, sel_factor=mzml_pheno_samples_pos, tune_components=length(unique(mzml_pheno_samples_pos)) - 1, sel_components=c(2,3), folds_number=8, keepx=model_splsda_keepx_pos, plot_roc_filename="plots/pos_ms2_classes_select_splsda_roc.pdf")
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_splsda_pos$'_selected_variables_')))
#f.heatmap.selected_features(feat_list=class_list_pos, sel_feat=sel_splsda_pos$'_selected_variables_', filename="plots/pos_ms2_classes_select_splsda.pdf", main="sPLS-DA")



# ############################## MS2 Superclasses statistics ##############################



# ---------- Histogram ----------
pdf(file="plots/pos_ms2_superclasses_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(unlist(superclass_list_pos)), main="Histogram", xlab="superclass_list")
dev.off()



# ---------- Variation partitioning ----------
model_varpart_pos <- varpart(superclass_list_pos, ~ mzml_pheno_samples_pos, ~ mzml_pheno_samples_pos)

# Plot results
pdf(file="plots/pos_ms2_superclasses_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart_pos, Xnames=c("sample group","samples"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca_pos <- prcomp(superclass_list_pos)

pdf(paste("plots/pos_ms2_superclasses_pca.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca_pos$x[, 2], model_pca_pos$x[,3], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca_pos)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, bg=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(model_pca_pos$x[, 2], model_pca_pos$x[,3], labels=mzml_names_pos, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()



# ---------- Feature selection of superclasses ----------
# Tukey on ANOVA
#sel_aov_pos <- f.select_features_aov(feat_class=mzml_pheno_samples_pos, feat_list=superclass_list_pos, conf_level=0.95, pvalue=0.99)
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_aov_pos)))
#f.heatmap.selected_features(feat_list=superclass_list_pos, sel_feat=sel_aov_pos, filename="plots/pos_ms2_superclasses_select_aov.pdf", main="Tukey on ANOVA")

# Linear Model
#sel_lm_pos <- f.select_features_linear_model(feat_matrix=superclass_list_pos, sel_factor=mzml_pheno_samples_pos, sel_formula=feature ~ mzml_pheno_samples_pos, pvalue=0.05)
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_lm_pos)))
#f.heatmap.selected_features(feat_list=superclass_list_pos, sel_feat=sel_lm_pos, filename="plots/pos_ms2_superclasses_select_lm.pdf", main="Linear Model")

# Random Forest
sel_rf_pos <- f.select_features_random_forest(feat_matrix=superclass_list_pos, sel_factor=mzml_pheno_samples_pos, sel_colors=mzml_pheno_colors_pos$cols, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/pos_ms2_superclasses_select_rf_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_rf_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=superclass_list_pos, sel_feat=sel_rf_pos$`_selected_variables_`, filename="plots/pos_ms2_superclasses_select_rf.pdf", main="Random Forest")

# PLS
sel_pls_pos <- f.select_features_pls(feat_matrix=superclass_list_pos, sel_factor=mzml_pheno_samples_pos, sel_colors=mzml_pheno_colors_pos$cols, components=2, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/pos_ms2_superclasses_select_pls_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_pls_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=superclass_list_pos, sel_feat=sel_pls_pos$`_selected_variables_`, filename="plots/pos_ms2_superclasses_select_pls.pdf", main="PLS")

# sPLS-DA
# First: Variation partitioning
#model_varpart_pos <- varpart(scale(superclass_list_pos), ~ mzml_pheno_samples_pos, ~ mzml_pheno_samples_pos)

# 10% of features correlate with factor
#model_varpart_corr_pos <- trunc(model_varpart_pos$part$indfract$Adj.R.squared[2] * ncol(superclass_list_pos))
#model_splsda_keepx_pos <- trunc(seq(model_varpart_corr_pos, 1,length.out=length(unique(mzml_pheno_samples_pos))))

#sel_splsda_pos <- f.select_features_splsda(feat_matrix=superclass_list_pos, sel_colors=mzml_pheno_colors_pos$cols, sel_factor=mzml_pheno_samples_pos, tune_components=length(unique(mzml_pheno_samples_pos)) - 1, sel_components=c(2,3), folds_number=8, keepx=model_splsda_keepx_pos, plot_roc_filename="plots/pos_ms2_superclasses_select_splsda_roc.pdf")
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_splsda_pos$'_selected_variables_')))
#f.heatmap.selected_features(feat_list=superclass_list_pos, sel_feat=sel_splsda_pos$'_selected_variables_', filename="plots/pos_ms2_superclasses_select_splsda.pdf", main="sPLS-DA")


