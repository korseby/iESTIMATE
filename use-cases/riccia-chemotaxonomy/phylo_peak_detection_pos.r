


# ############################## MS1 ##############################



# ---------- MS1 Preparations ----------
# MS1 variables
polarity <- "positive"
pol <- substr(x=polarity, start=1, stop=3)
ppm <- 35
phylo_ms1_intensity_cutoff <- 6500	          #approx. 0.01%

# General variables
phylo_mzml_files_pos <- NULL
phylo_mzml_names_pos <- NULL
phylo_mzml_times_pos <- NULL

# Load files
phylo_mzml_files_pos <- list.files(mzml_dir, pattern="*.mzML", recursive=T, full.names=T)
phylo_mzml_files_pos <- phylo_mzml_files_pos[grep(pol, phylo_mzml_files_pos, invert=FALSE)]
phylo_mzml_files_pos <- phylo_mzml_files_pos[grep("(R-|Lunularia)", phylo_mzml_files_pos, invert=FALSE)]
# Unfortunately, R-wanstorfii-1 is mis-identification of R-glauca !
phylo_mzml_files_pos <- phylo_mzml_files_pos[grep("R-warnstorfii-1", phylo_mzml_files_pos, invert=TRUE)]
phylo_mzml_files_pos <- phylo_mzml_files_pos[grep("R-bifurca", phylo_mzml_files_pos, invert=TRUE)]

# Basenames of files without path and without extension
phylo_mzml_names_pos <- gsub('(.*)\\..*', '\\1', gsub('( |-|,)', '.', basename(phylo_mzml_files_pos)))

# Create phenodata based on species
phylo_mzml_pheno_pos <- data.frame(sample_name=phylo_mzml_names_pos, sample_group=gsub(pattern="\\.\\d.*", replacement="", x=phylo_mzml_names_pos, perl=TRUE))
phylo_mzml_pheno_samples_pos <- as.factor(phylo_mzml_pheno_pos$sample_group)
phylo_mzml_pheno_colors_pos <- c("mediumorchid4","goldenrod4","skyblue4","forestgreen")
phylo_mzml_pheno_colors_samples_pos <- c("mediumorchid4","mediumorchid4","mediumorchid4","goldenrod4","goldenrod4","goldenrod4","skyblue4","skyblue4","skyblue4","forestgreen","forestgreen","forestgreen")

# Save timestamps of samples
for (i in 1:length(phylo_mzml_files_pos)) {
	fl <- mzR::openMSfile(phylo_mzml_files_pos[i])
	run_info <- mzR::runInfo(fl)
	mzR::close(fl)
	phylo_mzml_times_pos <- c(phylo_mzml_times_pos, run_info$startTimeStamp)
}

# Display MSn levels
phylo_mzml_msn_pos <- NULL
for (i in 1:length(phylo_mzml_files_pos)) {
	phylo_mzml_data_pos <- readMSData(phylo_mzml_files_pos[i], mode="onDisk")
	phylo_mzml_msn_pos <- rbind(phylo_mzml_msn_pos, t(as.matrix(table(msLevel(phylo_mzml_data_pos)))))
}
colnames(phylo_mzml_msn_pos) <- c("MS1", "MS2")
rownames(phylo_mzml_msn_pos) <- phylo_mzml_names_pos

# Plot MSn levels
pdf(file="plots/phylo_pos_msn_levels.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=16, family="Helvetica")
par(mfrow=c(2,1), mar=c(16,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
boxplot(phylo_mzml_msn_pos, main="Number of spectra")

model_boxplot <- boxplot(t(phylo_mzml_msn_pos[,2]), main="Number of MS2 spectra per sample", xaxt="n")
tick <- seq_along(model_boxplot$names)
axis(1, at=tick, labels=F)
text(tick, par("usr")[3]-par("usr")[3]/10, model_boxplot$names, adj=0, srt=270, xpd=T)
dev.off()



# ---------- Peak detection ----------
# Import raw data as MSnbase object
phylo_raw_data_pos <- readMSData(files=phylo_mzml_files_pos, pdata=new("NAnnotatedDataFrame", phylo_mzml_pheno_pos), mode="onDisk", centroided=TRUE)
table(msLevel(phylo_raw_data_pos))
head(fData(phylo_raw_data_pos)[, c("isolationWindowTargetMZ", "isolationWindowLowerOffset",
								   "isolationWindowUpperOffset", "msLevel", "retentionTime")])
write.csv(fData(phylo_raw_data_pos), file="data/phylo_pos_raw_data.csv", row.names=FALSE)

# Restrict data to 1020 seconds (17 minutes)
phylo_raw_data_pos <- filterRt(phylo_raw_data_pos, c(0, 1020))

# Inspect mz values per file
phylo_raw_mz_pos <- mz(phylo_raw_data_pos)
phylo_raw_mz_pos <- split(phylo_raw_mz_pos, f = fromFile(phylo_raw_data_pos))
print(length(phylo_raw_mz_pos))

# Get base peak chromatograms
phylo_chromas_pos <- chromatogram(phylo_raw_data_pos, aggregationFun="max")

# Plot chromatograms based on phenodata groups
pdf(file="plots/phylo_pos_chromas.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(phylo_chromas_pos, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col=phylo_mzml_pheno_colors_samples_pos)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=phylo_mzml_pheno_colors_pos, legend=phylo_mzml_pheno_samples_pos)
dev.off()

# Get TICs
pdf(file="plots/phylo_pos_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(5,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
tics_pos <- split(tic(phylo_raw_data_pos), f=fromFile(phylo_raw_data_pos))
boxplot(tics_pos, names=phylo_mzml_names_pos, las=2, col=phylo_mzml_pheno_colors_pos, ylab="intensity", main="Total ion current")
dev.off()

# Grouping/binning the samples based on similarity of their base peak chromatogram to spot potentially problematic samples
phylo_chromas_bin_pos <- MSnbase::bin(phylo_chromas_pos, binSize=2)
phylo_chromas_bin_cor_pos <- cor(log2(do.call(cbind, lapply(phylo_chromas_bin_pos, intensity))))
colnames(phylo_chromas_bin_cor_pos) <- rownames(phylo_chromas_bin_cor_pos) <- phylo_raw_data_pos$sample_name
pdf(file="plots/phylo_pos_chromas_bin_cor.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(phylo_chromas_bin_cor_pos)
dev.off()

# Assess retention times and intensities of first file
head(rtime(phylo_chromas_pos[1, 1]))
head(intensity(phylo_chromas_pos[1, 1]))

# Peak detection in MS1 data
if (polarity=="negative") {
	phylo_ms1_params_pos <- CentWaveParam(ppm=9.5, mzCenterFun="mean", peakwidth=c(4, 36), prefilter=c(2, 170), mzdiff=0.0045, snthresh=11, noise=0, integrate=1,
										  firstBaselineCheck=TRUE, verboseColumns=TRUE, fitgauss=FALSE, roiList=list(), roiScales=numeric())
} else {
	phylo_ms1_params_pos <- CentWaveParam(ppm=9.5, mzCenterFun="mean", peakwidth=c(4, 21), prefilter=c(2, 100), mzdiff=0.0034, snthresh=11, noise=0, integrate=1,
										  firstBaselineCheck=TRUE, verboseColumns=FALSE, fitgauss=FALSE, roiList=list(), roiScales=numeric())
}
phylo_ms1_data_pos <- findChromPeaks(phylo_raw_data_pos, param=phylo_ms1_params_pos)

# Per file summary
phylo_ms1_summary_pos <- lapply(split.data.frame(chromPeaks(phylo_ms1_data_pos), f=chromPeaks(phylo_ms1_data_pos)[, "sample"]), FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )
phylo_ms1_summary_pos <- do.call(rbind, phylo_ms1_summary_pos)
rownames(phylo_ms1_summary_pos) <- basename(fileNames(phylo_ms1_data_pos))
print(phylo_ms1_summary_pos)
table(msLevel(phylo_ms1_data_pos))
write.csv(as.data.frame(table(msLevel(phylo_ms1_data_pos))), file="data/phylo_pos_ms1_data.csv", row.names=FALSE)

# To get a global overview of the peak detection we can plot the frequency of identified peaks per file along the retention time axis. This allows to identify time periods along the MS run in which a higher number of peaks was identified and evaluate whether this is consistent across files.
pdf(file="plots/phylo_pos_ms1_data.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,18,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plotChromPeakImage(phylo_ms1_data_pos, main="Frequency of identified peaks per RT")
dev.off()

# Group peaks
if (polarity=="negative") {
	phylo_ms1_data_pos <- groupChromPeaks(phylo_ms1_data_pos, param=PeakDensityParam(sampleGroups=phylo_ms1_data_pos$sample_group, minFraction=0.7, bw=0.25))
} else {
	phylo_ms1_data_pos <- groupChromPeaks(phylo_ms1_data_pos, param=PeakDensityParam(sampleGroups=phylo_ms1_data_pos$sample_group, minFraction=0.7, bw=0.25))
}

# RT correction
if (polarity=="negative") {
	phylo_ms1_data_pos <- adjustRtime(phylo_ms1_data_pos, param=PeakGroupsParam(minFraction=0.7,smooth="loess",span=0.2,family="gaussian"))
} else {
	phylo_ms1_data_pos <- adjustRtime(phylo_ms1_data_pos, param=PeakGroupsParam(minFraction=0.7,smooth="loess",span=0.2,family="gaussian"))
}

# Plot the difference of raw and adjusted retention times
pdf(file="plots/phylo_pos_ms1_phylo_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(phylo_chromas_pos, peakType="none", main="Raw chromatograms", col=phylo_mzml_pheno_colors_samples_pos)
plotAdjustedRtime(phylo_ms1_data_pos, lwd=2, main="Retention Time correction", col=phylo_mzml_pheno_colors_samples_pos)
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
dev.off()

# Group peaks
if (polarity=="negative") {
	phylo_ms1_data_pos <- groupChromPeaks(phylo_ms1_data_pos, param=PeakDensityParam(sampleGroups=phylo_ms1_data_pos$sample_group, minFraction=0.7, bw=0.25))
} else {
	phylo_ms1_data_pos <- groupChromPeaks(phylo_ms1_data_pos, param=PeakDensityParam(sampleGroups=phylo_ms1_data_pos$sample_group, minFraction=0.7, bw=0.25))
}

# Get integrated peak intensity per feature/sample
print(head(featureValues(phylo_ms1_data_pos, value="into")))

# Fill peaks
#phylo_ms1_data_pos <- fillChromPeaks(phylo_ms1_data_pos, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
#head(featureValues(phylo_ms1_data_pos))
#head(featureSummary(phylo_ms1_data_pos, group=phylo_ms1_data_pos$sample_group))

# Evaluate grouping
pdf(file="plots/phylo_pos_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
phylo_ms1_pca_pos <- prcomp(t(na.omit(log2(featureValues(phylo_ms1_data_pos, value="into")))), center=TRUE)
plot(phylo_ms1_pca_pos$x[, 2], phylo_ms1_pca_pos$x[,3], pch=19, main="PCA: Grouping of samples",
	 xlab=paste0("PC1: ", format(summary(phylo_ms1_pca_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(phylo_ms1_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples_pos, cex=2)
grid()
text(phylo_ms1_pca_pos$x[, 2], phylo_ms1_pca_pos$x[,3], labels=phylo_ms1_data_pos$sample_name, col=phylo_mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()

# Show peaks
tail(chromPeaks(phylo_ms1_data_pos))
tail(chromPeakData(phylo_ms1_data_pos))

# Show process history
processHistory(phylo_ms1_data_pos)



# ---------- Build MS1 feature tables ----------
# Build feature matrix
phylo_ms1_matrix_pos <- featureValues(phylo_ms1_data_pos, method="medret", value="into")
colnames(phylo_ms1_matrix_pos) <- phylo_mzml_names_pos
dim(phylo_ms1_matrix_pos)
phylo_feat_list_pos <- t(phylo_ms1_matrix_pos)

# Build feature summary
phylo_ms1_summary_pos <- featureSummary(phylo_ms1_data_pos)
phylo_ms1_def_pos <- featureDefinitions(phylo_ms1_data_pos)

# Missing value imputation
#phylo_feat_list_pos[is.na(phylo_feat_list_pos)] <- median(na.omit(as.numeric(unlist(phylo_feat_list_pos))))

# Transform data
phylo_feat_list_pos <- log2(phylo_feat_list_pos)

# Missing value imputation
#phylo_feat_list_pos[which(is.na(phylo_feat_list_pos))] <- median(na.omit(as.numeric(unlist(phylo_feat_list_pos))))

phylo_feat_list_pos[which(is.na(phylo_feat_list_pos))] <- runif(length(which(is.na(phylo_feat_list_pos))), min=0, max=0.0000001)
#phylo_feat_list_pos <- scale(phylo_feat_list_pos, scale=FALSE, center=FALSE)
#phylo_feat_list_pos[is.na(phylo_feat_list_pos)] <- 0
#phylo_feat_list_pos[which(phylo_feat_list_pos < 0)] <- 0
#phylo_feat_list_pos[is.infinite(phylo_feat_list_pos)] <- 0
#phylo_feat_list_pos <- phylo_feat_list_pos[!apply(phylo_feat_list_pos, MARGIN=1, function(x) max(x,na.rm=TRUE) == min(x,na.rm=TRUE)),]

# Plot histogram
pdf(file="plots/phylo_pos_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(phylo_feat_list_pos), main="Histogram of feature table")
dev.off()

# PCA
pdf(file="plots/phylo_pos_ms1_feature_table_pca12.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
phylo_ms1_pca_pos <- prcomp(phylo_feat_list_pos, center=TRUE)
plot(phylo_ms1_pca_pos$x[, 1], phylo_ms1_pca_pos$x[,2], pch=19, main="PCA of feature table",
	 xlab=paste0("PC1: ", format(summary(phylo_ms1_pca_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(phylo_ms1_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples_pos, cex=2)
grid()
text(phylo_ms1_pca_pos$x[, 1], phylo_ms1_pca_pos$x[, 2], labels=phylo_ms1_data_pos$sample_name, col=phylo_mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()

pdf(file="plots/phylo_pos_ms1_feature_table_pca23.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
phylo_ms1_pca_pos <- prcomp(phylo_feat_list_pos, center=TRUE)
plot(phylo_ms1_pca_pos$x[, 2], phylo_ms1_pca_pos$x[,3], pch=19, main="PCA of feature table",
	 xlab=paste0("PC3: ", format(summary(phylo_ms1_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC4: ", format(summary(phylo_ms1_pca_pos)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=phylo_mzml_pheno_colors_samples_pos, cex=2)
grid()
text(phylo_ms1_pca_pos$x[, 2], phylo_ms1_pca_pos$x[, 3], labels=phylo_ms1_data_pos$sample_name, col=phylo_mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()

# Create single 0/1 matrix
phylo_bina_list_pos <- t(phylo_ms1_matrix_pos)
phylo_bina_list_pos[is.na(phylo_bina_list_pos)] <- 1
phylo_bina_list_pos <- log2(phylo_bina_list_pos)
phylo_bina_list_pos[phylo_bina_list_pos < log2(phylo_ms1_intensity_cutoff)] <- 0
phylo_bina_list_pos[phylo_bina_list_pos != 0] <- 1

# Only unique compounds in group phylo_mzml_pheno$ and not the others
phylo_uniq_list_pos <- apply(X=phylo_bina_list_pos, MARGIN=2, FUN=function(x) { if (length(unique(phylo_mzml_pheno_pos$sample_group[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
colnames(phylo_uniq_list_pos) <- colnames(phylo_bina_list_pos)
rownames(phylo_uniq_list_pos) <- rownames(phylo_bina_list_pos)

# Create data frame
phylo_model_div_pos             <- data.frame(features=apply(X=phylo_bina_list_pos, MARGIN=1, FUN=function(x) { sum(x) } ))
phylo_model_div_pos$richness    <- apply(X=phylo_bina_list_pos, MARGIN=1, FUN=function(x) { sum(x) } )
#phylo_model_div_pos$menhinick   <- apply(X=phylo_bina_list_pos, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
phylo_model_div_pos$shannon     <- apply(X=phylo_feat_list_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
phylo_model_div_pos$pielou      <- apply(X=scale(phylo_feat_list_pos, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#phylo_model_div_pos$chao        <- vegan::specpool2vect(X=vegan::specpool(phylo_feat_list_pos, species), index="chao")
phylo_model_div_pos$simpson     <- apply(X=phylo_feat_list_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
phylo_model_div_pos$inverse     <- apply(X=phylo_feat_list_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
phylo_model_div_pos$fisher      <- apply(X=phylo_feat_list_pos, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
phylo_model_div_pos$unique      <- apply(X=phylo_uniq_list_pos, MARGIN=1, FUN=function(x) { sum(x) })

# Remove NAs if present
phylo_model_div_pos[is.na(phylo_model_div_pos)] <- 0



# ############################## MS2 ##############################



# ---------- MS2 spectra detection ----------
# Estimate precursor intensity
phylo_precursor_intensity_pos <- estimatePrecursorIntensity(phylo_ms1_data_pos)
print(head(na.omit(phylo_precursor_intensity_pos)))

# Reconstruct MS2 spectra from MS1 data
phylo_ms2_data_pos <- chromPeakSpectra(phylo_ms1_data_pos, msLevel=2L, return.type="Spectra")
print(phylo_ms2_data_pos)
print(length(phylo_ms2_data_pos$peak_id))

# Extract all usable MS2 spectra
phylo_ms2_spectra_pos <- list()
#for (i in 1:nrow(phylo_ms1_def_pos)) {
phylo_ms2_spectra_pos <- foreach(i=1:nrow(phylo_ms1_def_pos)) %dopar% {
	# Extract existing MS2 spectra for feature
	feature_of_interest <- phylo_ms1_def_pos[i, "mzmed"]
	peaks_of_interest <- chromPeaks(phylo_ms1_data_pos, mz=feature_of_interest, ppm=ppm)
	
	# Continue if feature has MS2 peaks
	if (length(which(phylo_ms2_data_pos$peak_id %in% rownames(peaks_of_interest))) > 0) {
		# Extract spectra
		spectra_of_interest <- phylo_ms2_data_pos[phylo_ms2_data_pos$peak_id %in% rownames(peaks_of_interest)]
		combined_spectra_of_interest <- filterIntensity(spectra_of_interest, intensity=c(1,Inf), backend=MsBackendDataFrame)
		combined_spectra_of_interest <- setBackend(combined_spectra_of_interest, backend=MsBackendDataFrame())
		
		# Combine spectra
		combined_spectra_of_interest <- Spectra::combineSpectra(combined_spectra_of_interest, FUN=combinePeaks, ppm=ppm, peaks="union", minProp=0.8, intensityFun=median, mzFun=median, backend=MsBackendDataFrame)#f=rownames(peaks_of_interest))
		
		# Remove noise from spectra
		#combined_spectra_of_interest <- pickPeaks(combined_spectra_of_interest, snr=1.0, method="SuperSmoother") #MAD
		#combined_spectra_of_interest <- Spectra::smooth(combined_spectra_of_interest, method="SavitzkyGolay") #(Weighted)MovingAverage
		
		# Only keep spectral data
		combined_spectra_peaks <- as.data.frame(Spectra::peaksData(combined_spectra_of_interest)[[1]])
		
		# Plot merged spectrum
		#Spectra::plotSpectra(combined_spectra_of_interest)
		#plot(x=combined_spectra_peaks[,1], y=combined_spectra_peaks[,2], type="h", xlab="m/z", ylab="intensity", main=paste("Precursor m/z",combined_spectra_of_interest@backend@spectraData$precursorMz[[1]]))
		#length(spectra_of_interest$peak_id)
		
		#phylo_ms2_spectra_pos[[i]] <- combined_spectra_of_interest
		return(combined_spectra_of_interest)
	}
}

# Remove empty spectra
names(phylo_ms2_spectra_pos) <- rownames(phylo_ms1_def_pos)
phylo_ms2_spectra_pos <- phylo_ms2_spectra_pos[lengths(phylo_ms2_spectra_pos) != 0]

# Relate all MS2 spectra to MS1 precursors
phylo_ms1_def_pos$has_ms2 <- as.integer(rownames(phylo_ms1_def_pos) %in% names(phylo_ms2_spectra_pos))
print(paste0("Number of MS2 spectra related to precursor: ", length(which(phylo_ms1_def_pos$has_ms2>0))))



# ---------- Annotate MS2 spectra ----------
# Save all MS2 spectra in MSP file
msp_text <- NULL
alignment_id <- 0
for (i in names(phylo_ms2_spectra_pos)) {
	alignment_id <- alignment_id + 1
	
	msp_text <- c(msp_text, paste("NAME:", i))
	msp_text <- c(msp_text, paste("AlignmentID:", alignment_id))
	msp_text <- c(msp_text, paste("RETENTIONTIME:", phylo_ms1_def_pos[i, "rtmed"]))
	msp_text <- c(msp_text, paste("PRECURSORMZ:", phylo_ms1_def_pos[i, "mzmed"]))
	msp_text <- c(msp_text, paste("METABOLITENAME:", i))
	if (polarity == "positive") {
		msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
	} else {
		msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
	}
	msp_text <- c(msp_text, paste("NumPeaks:", nrow(as.data.frame(peaksData(phylo_ms2_spectra_pos[[i]])[[1]]))))
	msp_text <- c(msp_text, paste(as.data.frame(peaksData(phylo_ms2_spectra_pos[[i]])[[1]])$mz, as.data.frame(peaksData(phylo_ms2_spectra_pos[[i]])[[1]])$intensity, sep="\t"))
	msp_text <- c(msp_text, "")
	
}

# Write MSP file
cat(msp_text, file="data/phylo_pos_ms2_spectra.msp", sep="\n")


# Save all MS2 spectra in MGF file
mgf_text <- NULL
for (i in names(phylo_ms2_spectra_pos)) {
	mgf_text <- c(mgf_text, paste0("COM=", i))
	mgf_text <- c(mgf_text, "BEGIN IONS")
	mgf_text <- c(mgf_text, "MSLEVEL=2")
	mgf_text <- c(mgf_text, paste0("TITLE=", i))
	mgf_text <- c(mgf_text, paste0("RTINSECONDS=", phylo_ms1_def_pos[i, "rtmed"]))
	mgf_text <- c(mgf_text, paste0("PEPMASS=", phylo_ms1_def_pos[i, "mzmed"]))
	if (polarity == "positive") {
		mgf_text <- c(mgf_text, paste0("CHARGE=", "1+"))
	} else {
		mgf_text <- c(mgf_text, paste0("CHARGE=", "1-"))
	}
	mgf_text <- c(mgf_text, paste(as.data.frame(peaksData(phylo_ms2_spectra_pos[[i]])[[1]])$mz, as.data.frame(peaksData(phylo_ms2_spectra_pos[[i]])[[1]])$intensity, sep=" "))
	mgf_text <- c(mgf_text, "END IONS")
	mgf_text <- c(mgf_text, "")
}

# Write MGF file
cat(mgf_text, file="data/phylo_pos_ms2_spectra.mgf", sep="\n")



# ---------- Classify MS2 spectra with CANOPUS ----------
# Apply classified classes onto feature table
phylo_ms1_def_pos$primary_class <- ""
phylo_ms1_def_pos$alternative_classes <- ""

# Read SIRIUS/CANOPUS classifier
if (SIRIUS_VERSION == 4) {
	phylo_classifier_canopus_pos <- read.table(file=paste0("data/pos_ms2_canopus_summary.tsv"), header=TRUE, sep="\t", quote="\"", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
	phylo_classifier_canopus_pos$Metabolite.name <- gsub(x=phylo_classifier_canopus_pos$name, pattern=".*_", replacement="")
	phylo_classifier_canopus_pos$primary_class <- paste("Organic compounds", phylo_classifier_canopus_pos$superclass, phylo_classifier_canopus_pos$class, phylo_classifier_canopus_pos$subclass, phylo_classifier_canopus_pos$level.5, sep="; ")
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	phylo_classifier_canopus_pos$Annotation..putative. <- phylo_classifier_canopus_pos$primary_class
	phylo_classifier_canopus_pos$alternative_classes <- phylo_classifier_canopus_pos$all.classifications
} else {
	phylo_classifier_canopus_pos <- read.table(file=paste0("data/pos_ms2_canopus_summary.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
	phylo_classifier_canopus_pos$Metabolite.name <- ""
	for (i in 1:length(phylo_classifier_canopus_pos$id)) {
		x = which(gsub(x=phylo_classifier_canopus_pos$id, pattern=".*?_", replacement="", perl=TRUE) %in% gsub(x=phylo_classifier_canopus_pos$id[i], pattern=".*?_", replacement="", perl=TRUE))
		if (length(x) > 0) phylo_classifier_canopus_pos$Metabolite.name[i] <- gsub(x=phylo_classifier_canopus_pos$id[x], pattern=".*?_", replacement="", perl=TRUE)
	}
	phylo_classifier_canopus_pos$Metabolite.name[phylo_classifier_canopus_pos$Metabolite.name == "null"] <- ""
	phylo_classifier_canopus_pos$primary_class <- paste("Organic compounds", phylo_classifier_canopus_pos$ClassyFire.superclass, phylo_classifier_canopus_pos$ClassyFire.class, phylo_classifier_canopus_pos$ClassyFire.subclass, phylo_classifier_canopus_pos$ClassyFire.level.5, sep="; ")
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	phylo_classifier_canopus_pos$Annotation..putative. <- phylo_classifier_canopus_pos$primary_class
	phylo_classifier_canopus_pos$alternative_classes <- phylo_classifier_canopus_pos$all.classifications
}

for (j in unique(phylo_classifier_canopus_pos$Metabolite.name)) {
	obj <- phylo_classifier_canopus_pos[phylo_classifier_canopus_pos$Metabolite.name %in% j, "Annotation..putative."]
	phylo_ms1_def_pos[j, "primary_class"] <- obj[1]
}



# ---------- Diversity of MS2 classes ----------
# Create CANOPUS classifier object for each sample
phylo_classifiers_pos <- list()
for (i in phylo_mzml_names_pos) {
	obj <- names(which(phylo_bina_list_pos[rownames(phylo_bina_list_pos)==i, colnames(phylo_bina_list_pos) %in% phylo_classifier_canopus_pos$Metabolite.name] > 0))
	phylo_classifiers_pos[[i]] <- phylo_classifier_canopus_pos[phylo_classifier_canopus_pos$Metabolite.name %in% obj, ]
}

# Diversity of classes per sample
phylo_div_classes_samples_pos <- NULL
for (i in phylo_mzml_names_pos) {
	obj <- table(phylo_classifiers_pos[[i]][,"Annotation..putative."])
	obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
	if (is.null(phylo_div_classes_samples_pos)) {
		phylo_div_classes_samples_pos <- obj
	} else {
		phylo_div_classes_samples_pos <- merge(phylo_div_classes_samples_pos, obj, by="classes", all.x=TRUE, all.y=TRUE)
	}
}
rownames(phylo_div_classes_samples_pos) <- phylo_div_classes_samples_pos$classes
phylo_div_classes_samples_pos <- phylo_div_classes_samples_pos[, -which(colnames(phylo_div_classes_samples_pos)=="classes")]
colnames(phylo_div_classes_samples_pos) <- phylo_mzml_names_pos

# Diversity of classes
phylo_div_classes_pos <- phylo_div_classes_samples_pos
phylo_div_classes_pos[is.na(phylo_div_classes_pos)] <- 0
phylo_div_classes_pos <- apply(X=phylo_div_classes_pos, MARGIN=1, FUN=function(x) { sum(x) })
phylo_div_classes_pos <- data.frame(row.names=names(phylo_div_classes_pos), frequency=as.numeric(phylo_div_classes_pos))

# Plot diversity of classes in all samples
pdf(file="plots/phylo_pos_ms2_classes_diversity.pdf", encoding="ISOLatin1", pointsize=8, width=6, height=14, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,15,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(phylo_div_classes_pos$frequency, names.arg=gsub('.*; ','',rownames(phylo_div_classes_pos)), las=1, horiz=TRUE, xlab="frequency", main="Diversity of compound classes", col=rainbow(n=nrow(phylo_div_classes_pos), alpha=0.6))
dev.off()

# Sunburst plot of classes of all samples
pdf(file="plots/phylo_pos_ms2_classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(phylo_div_classes_pos), phylo_div_classes_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(phylo_div_classes_pos), "Number of spectra"=phylo_div_classes_pos$frequency), file="plots/phylo_pos_ms2_classes_sunburst.csv", row.names=FALSE)

# Classes
phylo_classes_pos <- rownames(phylo_div_classes_pos)
phylo_classes_pos <- phylo_classes_pos[which(grepl(pattern="^Organic compounds", x=phylo_classes_pos))]
phylo_classes_pos <- gsub(x=phylo_classes_pos, pattern=":", replacement="")
phylo_classes_pos <- gsub(x=phylo_classes_pos, pattern="/", replacement="; ")



# ---------- Build diversity objects ----------
# Imputation of NA with zeros
phylo_div_classes_pos[is.na(phylo_div_classes_pos)] <- 0
phylo_div_classes_samples_pos[is.na(phylo_div_classes_samples_pos)] <- 0

# Classification list for statistics
phylo_class_list_pos <- as.data.frame(t(phylo_div_classes_samples_pos))
phylo_class_list_pos[is.na(phylo_class_list_pos)] <- 0

# Log Transform
#phylo_class_list_pos <- log2(phylo_class_list_pos + 1)

# Only keep class names
colnames(phylo_class_list_pos) <- gsub(x=colnames(phylo_class_list_pos), pattern='.*; ', replacement='')

# Generate phylo_class_int_list_pos with abundances instead of counts
phylo_class_int_list_pos <- phylo_class_list_pos

for (i in 1:nrow(phylo_class_list_pos)) {
	samp <- rownames(phylo_class_list_pos)[i]
	for (j in 1:ncol(phylo_class_list_pos)) {
		cl <- colnames(phylo_class_list_pos)[j]
		ft <- phylo_classifier_canopus_pos$Metabolite.name[which(gsub(x=phylo_classifier_canopus_pos$primary_class, pattern='.*; ', replacement='') == cl)]
		ints <- phylo_feat_list_pos[i, which(colnames(phylo_feat_list_pos) %in% ft)]
		phylo_class_int_list_pos[i, j] <- sum(ints) * phylo_class_list_pos[i, j]
	}
}



# ---------- Build diversity objects for in- and outgroup ----------
# Riccia
phylo_riccia_div_classes_samples_pos <- phylo_div_classes_samples_pos[,4:12]

phylo_riccia_div_classes_pos <- phylo_riccia_div_classes_samples_pos
phylo_riccia_div_classes_pos[is.na(phylo_riccia_div_classes_pos)] <- 0
phylo_riccia_div_classes_pos <- apply(X=phylo_riccia_div_classes_pos, MARGIN=1, FUN=function(x) { sum(x) })
phylo_riccia_div_classes_pos <- data.frame(row.names=names(phylo_riccia_div_classes_pos), frequency=as.numeric(phylo_riccia_div_classes_pos))

phylo_riccia_classes_pos <- rownames(phylo_riccia_div_classes_pos)
phylo_riccia_classes_pos <- phylo_riccia_classes_pos[which(grepl(pattern="^Organic compounds", x=phylo_riccia_classes_pos))]
phylo_riccia_classes_pos <- gsub(x=phylo_riccia_classes_pos, pattern=":", replacement="")
phylo_riccia_classes_pos <- gsub(x=phylo_riccia_classes_pos, pattern="/", replacement="; ")

# Lunularia
phylo_lunularia_div_classes_samples_pos <- phylo_div_classes_samples_pos[,1:3]

phylo_lunularia_div_classes_pos <- phylo_lunularia_div_classes_samples_pos
phylo_lunularia_div_classes_pos[is.na(phylo_lunularia_div_classes_pos)] <- 0
phylo_lunularia_div_classes_pos <- apply(X=phylo_lunularia_div_classes_pos, MARGIN=1, FUN=function(x) { sum(x) })
phylo_lunularia_div_classes_pos <- data.frame(row.names=names(phylo_lunularia_div_classes_pos), frequency=as.numeric(phylo_lunularia_div_classes_pos))

phylo_lunularia_classes_pos <- rownames(phylo_lunularia_div_classes_pos)
phylo_lunularia_classes_pos <- phylo_lunularia_classes_pos[which(grepl(pattern="^Organic compounds", x=phylo_lunularia_classes_pos))]
phylo_lunularia_classes_pos <- gsub(x=phylo_lunularia_classes_pos, pattern=":", replacement="")
phylo_lunularia_classes_pos <- gsub(x=phylo_lunularia_classes_pos, pattern="/", replacement="; ")



# ---------- Classification at level of superclasses ----------
if (FALSE) {
	# Use CANOPUS Metabolite.name classification
	phylo_classifier_canopus_pos$primary_class <- paste("Organic compounds", phylo_classifier_canopus_pos$superclass, phylo_classifier_canopus_pos$most.specific.class, sep="; ")
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	phylo_classifier_canopus_pos$primary_class <- gsub(x=phylo_classifier_canopus_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	phylo_classifier_canopus_pos$Annotation..putative. <- phylo_classifier_canopus_pos$primary_class
	
	# Create CANOPUS classifier object for each sample
	phylo_superclassifiers_pos <- list()
	for (i in phylo_mzml_names_pos) {
		obj <- names(which(phylo_bina_list_pos[rownames(phylo_bina_list_pos)==i, colnames(phylo_bina_list_pos) %in% phylo_classifier_canopus_pos$Metabolite.name] > 0))
		phylo_superclassifiers_pos[[i]] <- phylo_classifier_canopus_pos[phylo_classifier_canopus_pos$Metabolite.name %in% obj, ]
	}
	
	# Diversity of classes per sample
	phylo_div_superclasses_samples_pos <- NULL
	for (i in phylo_mzml_names_pos) {
		obj <- table(phylo_superclassifiers_pos[[i]][,"Annotation..putative."])
		obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
		if (is.null(phylo_div_superclasses_samples_pos)) {
			phylo_div_superclasses_samples_pos <- obj
		} else {
			phylo_div_superclasses_samples_pos <- merge(phylo_div_superclasses_samples_pos, obj, by="classes", all.x=TRUE, all.y=TRUE)
		}
	}
	rownames(phylo_div_superclasses_samples_pos) <- phylo_div_superclasses_samples_pos$classes
	phylo_div_superclasses_samples_pos <- phylo_div_superclasses_samples_pos[, -which(colnames(phylo_div_superclasses_samples_pos)=="classes")]
	colnames(phylo_div_superclasses_samples_pos) <- phylo_mzml_names_pos
} else {
	# Make superclasses at level 3
	phylo_superclass_level_pos <- 3
	phylo_div_superclasses_samples_names_pos <- NULL
	for (i in c(1:phylo_superclass_level_pos)) {
		phylo_div_superclasses_samples_names_pos <- c(phylo_div_superclasses_samples_names_pos, lapply(X=strsplit(rownames(phylo_div_classes_samples_pos), '; '), FUN=function(x) { gsub(x=paste(x[1:i],sep='',collapse='; '),pattern='; NA',replacement='') }))
	}
	phylo_div_superclasses_samples_pos <- data.frame()
	for (i in c(1:ncol(phylo_div_classes_samples_pos))) phylo_div_superclasses_samples_pos <- rbind(phylo_div_superclasses_samples_pos, rep(0, length(unique(phylo_div_superclasses_samples_names_pos))))
	phylo_div_superclasses_samples_pos <- t(phylo_div_superclasses_samples_pos)
	colnames(phylo_div_superclasses_samples_pos) <- colnames(phylo_div_classes_samples_pos)
	rownames(phylo_div_superclasses_samples_pos) <- unique(phylo_div_superclasses_samples_names_pos)
	for (i in rownames(phylo_div_superclasses_samples_pos)) {
		for (j in c(1:ncol(phylo_div_classes_samples_pos))) {
			phylo_div_superclasses_samples_pos[rownames(phylo_div_superclasses_samples_pos)==i, j] <- sum(phylo_div_classes_samples_pos[grep(x=rownames(phylo_div_classes_samples_pos), pattern=i), j])
		}
	}
}

# Diversity of superclasses
phylo_div_superclasses_pos <- phylo_div_superclasses_samples_pos
phylo_div_superclasses_pos[is.na(phylo_div_superclasses_pos)] <- 0
phylo_div_superclasses_pos <- apply(X=phylo_div_superclasses_pos, MARGIN=1, FUN=function(x) { sum(x) })
phylo_div_superclasses_pos <- data.frame(row.names=names(phylo_div_superclasses_pos), frequency=as.numeric(phylo_div_superclasses_pos))

# Sunburst plot
pdf(file="plots/phylo_pos_ms2_superclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(phylo_div_superclasses_pos), phylo_div_superclasses_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(phylo_div_superclasses_pos), "Number of spectra"=phylo_div_superclasses_pos$frequency), file="plots/phylo_pos_ms2_superclasses_sunburst.csv", row.names=FALSE)

# Imputation of NA with zeros
phylo_div_superclasses_pos[is.na(phylo_div_superclasses_pos)] <- 0
phylo_div_superclasses_samples_pos[is.na(phylo_div_superclasses_samples_pos)] <- 0

# Classification list for statistics
phylo_superclass_list_pos <- as.data.frame(t(phylo_div_superclasses_samples_pos))
phylo_superclass_list_pos[is.na(phylo_superclass_list_pos)] <- 0

# Log Transform
#phylo_superclass_list_pos <- log2(phylo_superclass_list_pos + 1)

# Only keep superclass names
colnames(phylo_superclass_list_pos) <- gsub(x=colnames(phylo_superclass_list_pos), pattern='.*; ', replacement='')

# Generate phylo_superclass_int_list_pos with abundances instead of counts
phylo_superclass_int_list_pos <- phylo_superclass_list_pos

for (i in 1:nrow(phylo_superclass_list_pos)) {
	samp <- rownames(phylo_superclass_list_pos)[i]
	for (j in 1:ncol(phylo_superclass_list_pos)) {
		cl <- colnames(phylo_superclass_list_pos)[j]
		ft <- phylo_classifier_canopus_pos$Metabolite.name[which(phylo_classifier_canopus_pos$primary_class %in% phylo_classifier_canopus_pos$primary_class[grep(x=phylo_classifier_canopus_pos$primary_class, pattern=cl)]) ]
		ints <- as.numeric(phylo_feat_list_pos[i, which(colnames(phylo_feat_list_pos) %in% ft)])
		phylo_superclass_int_list_pos[i, j] <- sum(ints) * as.numeric(phylo_superclass_list_pos[i, j])
	}
}


