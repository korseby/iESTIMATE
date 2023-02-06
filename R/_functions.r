
# ---------- Shannon Diversity ----------
shannon.diversity <- function(p) {
	# Based on Li et al. (2016)
	# Function is obsolete, as it returns same values than vegan::diversity(x, index="shannon")
	# Since we do not expect features with negative intensity,
	# we exclude negative values and take the natural logarithm
	if (min(p) < 0 || sum(p) <= 0) 
		return(NA)
	pij <- p[p>0] / sum(p) 
	-sum(pij * log(pij)) 
}



# ---------- Menhinick Diversity ----------
menhinick.diversity <- function(p) {
	# Based on: http://www.coastalwiki.org/wiki/Measurements_of_biodiversity#Species_richness_indices
	D_Mn <- length(p) / sqrt(vegan::specnumber(p))
}



# ---------- Tukey-Test ----------
tukey.test <- function(response, term) {
	model_anova <- aov(formula(response ~ term))
	model_mc <- multcomp::glht(model_anova, multcomp::mcp(term="Tukey"))
	model_cld <- multcomp::cld(summary(model_mc), decreasing=TRUE)
	model_tukey <- data.frame("tukey_groups"=model_cld$mcletters$Letters)
	return(model_tukey)
}



# ---------- P-Values ----------
print_p.values <- function(p.values) {
	p.values[1] <- 1
	p.values[p.values < 0.001] <- '***'
	p.values[(p.values >= 0.001 & p.values < 0.01)] <- '**'
	p.values[(p.values >= 0.01 & p.values < 0.05)] <- '*'
	p.values[(p.values >= 0.05 & p.values < 0.1)] <- '.'
	p.values[p.values >= 0.1] <- ' '
	return(p.values)
}



# ---------- Imputation using NIPALS ----------
f.impute_nipals <- function(feat_matrix, components) {
	# Impute NAs by performing NIPALS
	model_nipals <- mixOmics::nipals(X=feat_matrix, ncomp=components, reconst=TRUE, max.iter=1000, tol=1e-09)
	
	# Return imputed data matrix
	return(model_nipals$rec)
}



# ---------- Calculate R-squared ----------
f.r2 <- function(actual, predicted) {
	actual <- as.numeric(actual)
	predicted <- as.numeric(predicted)

	R2 <- caret::postResample(actual, predicted)[2]	
	#R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))
	
	return(R2)
}



# ---------- Calculate Classification Rate ----------
f.classification_rate <- function(sel_factor, predicted) {
	accuracy <- table(predicted, sel_factor)
	cr <- sum(diag(accuracy))/sum(accuracy)
	
	#print(paste0("Classification rate: ", round(cr ,3)))
	
	return(cr)
}



# ---------- Calculate Accuracy ----------
f.accuracy <- function(sel_factor, predicted) {
	return(length(which(predicted == sel_factor)) / length(sel_factor))
}



# ---------- Calculate Weighted Accuracy ----------
f.weighted_accuracy <- function(sel_factor, predicted) {
	# Calculate weights
	weights <- rep(1 / length(levels(sel_factor)), length(levels(sel_factor)))
	
	# Check
	sel_levels <- levels(sel_factor)
	if (length(weights) != length(sel_levels)) {
		stop("Error! Number of weights should have some length as the number of classes.")
	}
	if (sum(weights) != 1) {
		stop("Error! Weights do not sum to 1.")
	}
	
	# Calculate
	accuracies <- lapply(sel_levels, function(x) {
		idx <- which(sel_factor==x)
		return(f.accuracy(predicted[idx], sel_factor[idx]))
	})
	accuracies <- mean(unlist(accuracies))
	
	return(accuracies)
}



# ---------- Performance Measures ----------
f.performance_measures_caret <- function(model, sel_factor, sel_colors) {
	# List of measures
	sel <- list()
	
	# Select index with highest accuracy from prediction model, only when not using savePredictions="final"
	#sel_ind_acc <- which(model$results$Accuracy==max(model$results$Accuracy))
	#if (! is.null(names(model$pred$mtry))) {
	#	sel_ind <- model$pred$mtry == unique(model$pred$mtry)[sel_ind_acc]
	#} else {
	#	sel_ind <- c(1:nrow(model$pred))
	#}
	
	# Performance measures
	sel_levels <- levels(as.factor(sel_factor))
	
	# Build prediction matrix
	model_pred <- as.data.frame(model$pred[, which(colnames(model$pred) %in% sel_levels)])
	colnames(model_pred) <- paste0(colnames(model_pred), "_pred_")
	#model_pred <- scales::rescale(x=as.matrix(model_pred), to=c(0,1))
	pred_true <- data.frame(dummies::dummy(model$pred$obs))
	colnames(pred_true) <- paste0(sel_levels, "_true")
	model_pred <- cbind(pred_true, model_pred)
	
	# Create classification objects
	sel_obs <- model$pred$obs
	sel_pred <- model$pred$pred
	sel_prob <- as.numeric(as.factor(sel_obs)); for (i in sel_levels) sel_prob[which(sel_obs==i)] <- as.numeric(model$pred[which(model$pred$obs==i), which(colnames(model$pred)==i)])
	
	# Multi-class AUC as defined by Hand & Till (2001)
	if (nlevels(sel_factor) > 2) {
		multiclass_roc <- pROC::multiclass.roc(response=sel_obs, predictor=sel_prob, levels=sel_levels, percent=FALSE, print.auc=TRUE)
		multiclass_auc <- round(as.numeric(multiclass_roc$auc),3)
	}
	
	# Multi-class classification metrics
	sel_pred <- as.factor(sel_pred)
	#levels(sel_pred) <- levels(sel_obs)
	model_metrics <- mltest::ml_test(predicted=sel_pred, true=sel_obs)
	
	# Multi-class Classification rate (= 1 - error rate)
	mcr <- f.classification_rate(sel_factor=sel_obs, predicted=sel_pred)
	
	# Balanced Error Rate (BER) (= 1 - multi-class classification rate)
	ber <- mlr::measureBER(truth=sel_obs, response=sel_pred)
	
	# Accuracy
	accuracy <- f.accuracy(sel_factor=sel_obs, predicted=sel_pred)
	weighted_accuracy <- f.weighted_accuracy(sel_factor=as.factor(sel_obs), predicted=as.factor(sel_pred))
	
	# R2
	R2 <- f.r2(actual=as.numeric(as.factor(sel_obs)), predicted=as.numeric(as.factor(sel_pred)))
	
	# ROC curve (one vs. rest)
	model_roc <- plot.roc(as.numeric(as.factor(sel_obs)), as.numeric(as.factor(sel_pred)), main="Confidence intervals", percent=TRUE, ci=TRUE, print.auc=TRUE)
	model_roc_confidence <- ci.se(model_roc, specificities=seq(0, 100, 1), of="thresholds", thresholds="local maximas") # Confidence
	model_roc_sensitivity <- ci.se(model_roc, specificities=seq(0, 100, 10), of="thresholds", thresholds="local maximas") # Sensitivity
	model_roc_specificity <- ci.sp(model_roc, sensitivities=seq(0, 100, 10), of="thresholds", thresholds="local maximas") # Specificity
	
	# PR curve
	#model_pr <- pr.curve(scores.class0=sel_prob[sel_obs == sel_pred], scores.class1=sel_prob[sel_obs != sel_pred], curve=TRUE)
	
	# Save metrics
	sel[["_model_roc_"]] <- model_roc
	sel[["_model_roc_confidence_"]] <- model_roc_confidence
	sel[["_model_roc_sensitivity_"]] <- model_roc_sensitivity
	sel[["_model_roc_specificity_"]] <- model_roc_specificity
	sel[["_multiclass_metrics_"]] <- model_metrics
	if (nlevels(sel_factor) > 2) { sel[["_model_multiclass_roc_"]] <- multiclass_roc }
	if (nlevels(sel_factor) > 2) { sel[["_multiclass_auc_"]] <- multiclass_auc }
	sel[["_multi_class_rate_"]] <- mcr
	sel[["_balanced_error_rate_"]] <- ber
	sel[["_accuracy_"]] <- accuracy
	sel[["_weighted_accuracy_"]] <- weighted_accuracy
	sel[["_model_r2_"]] <- R2
	sel[["_AUC_"]] <- round(as.numeric(model_roc$auc)/100, 3)
	#sel[["_AUCPR_"]] <- round(model_pr$auc.integral, 3)
	
	# Calculate ROC and PR with multiROC
	model_roc <- multi_roc(model_pred, force_diag=TRUE)
	model_pr <- multi_pr(model_pred, force_diag=TRUE)
	
	# Plot ROC using multiROC
	plot_roc <- plot_roc_data(model_roc)
	plot(x=c(0,1), y=c(0,1), type="l", xlim=c(0,1), ylim=c(0,1), xlab="1 - Specificity", ylab="Sensitivity", main="ROC")
	for (i in sel_levels) {
		lines(x=1-plot_roc[which(plot_roc$Group==i),"Specificity"], y=plot_roc[which(plot_roc$Group==i),"Sensitivity"], col=sel_colors[which(sel_levels==i)], lwd=2)
	}
	legend("bottomright", legend=paste0(model_roc$Groups, " (AUC=", round(unlist(model_roc$AUC)[model_roc$Groups]*100, 2), "%)"), col=sel_colors, lwd=2, cex=0.75)
	
	# Plot ROC using pROC
	model_proc <- list()
	for (i in c(1:length(sel_levels))) {
		if (i==1) add=FALSE else add=TRUE
		model_proc[[as.character(sel_levels[i])]] <- plot.roc(model_pred[,i], model_pred[,i+length(sel_levels)], main="ROC curves of levels", col=sel_colors[i], percent=TRUE, ci=TRUE, print.auc=FALSE, add=add)
		plot(ci.se(model_proc[[as.character(sel_levels[i])]], specificities=seq(0, 100, 1), of="thresholds", thresholds="local maximas"), type="shape", col=rgb(col2rgb(sel_colors[i])[1]/255, col2rgb(sel_colors[i])[2]/255, col2rgb(sel_colors[i])[3]/255, alpha=0.3), lty=0, no.roc=TRUE)
	}
	legend("bottomright", legend=paste0(sel_levels, ", AUC: ", round(as.numeric(unlist(lapply(model_proc, function(x) x$auc))), 1), "%"), col=sel_colors, lwd=2, cex=0.75)
	
	# Plot Precision and Recall Curve using multiROC
	plot_pr <- plot_pr_data(model_pr)
	plot(x=c(0,1), y=c(1,0), type="n", xlim=c(0,1), ylim=c(0,1), xlab="Recall", ylab="Precision", main="PR")
	for (i in sel_levels) {
		lines(x=plot_pr[which(plot_pr$Group==i),"Recall"], y=plot_pr[which(plot_pr$Group==i),"Precision"], col=sel_colors[which(sel_levels==i)], lwd=2)
	}
	legend("bottomleft", legend=paste0(model_pr$Groups, " (AUC-PR=", round(unlist(model_pr$AUC)[model_pr$Groups]*100, 2), "%)"), col=sel_colors, lwd=2, cex=0.75)
	
	# Plot Precision and Recall Curve using pROC
	model_ppr <- list()
	plot(xlim=c(0,1), ylim=c(0,1), x=NULL, y=NULL, xlab="Recall", ylab="Precision", main="Precision Recall Curves of levels")
	for (i in 1:length(sel_levels)) {
		model_ppr[[as.character(sel_levels[i])]] <- pr.curve(scores.class0=model_pred[as.logical(model_pred[,i]),i+length(sel_levels)], scores.class1=model_pred[!as.logical(model_pred[,i]),i+length(sel_levels)], curve=TRUE)
		lines(x=model_ppr[[as.character(sel_levels[i])]]$curve[,1], y=model_ppr[[as.character(sel_levels[i])]]$curve[,2], xlab="Recall",ylab="Precision", t="l", col=sel_colors[i], lwd=2)
	}
	legend("bottomleft", legend=paste0(sel_levels, ", AUC-PR: ", round(as.numeric(unlist(lapply(model_ppr, function(x) x$auc.integral*100))), 1), "%"), col=sel_colors, lwd=2, cex=0.75)
	
	# Save metrics
	sel[["_AUCPR_"]] <- model_ppr[[1]]$auc.integral
	sel[["_model_pr_"]] <- model_pr
	sel[["_model_factor_roc_"]] <- model_proc
	sel[["_model_factor_pr_"]] <- model_ppr
	
	return(sel)
}



# ---------- Linear Model ----------
f.select_features_linear_model <- function(feat_matrix, sel_factor, sel_formula, pvalue=0.05) {
	if (length(unique(sel_factor)) <= 2 ) {
		sel_lm <- NULL
	} else {
		print("Warning. Selection with more than 2 factors is not supported.")
		sel_lm <- list()
	}
	
	for (i in 1:ncol(feat_matrix)) {
		# Create data frame for feature of interest
		lm_list <- data.frame(feature=feat_matrix[,i])
		
		# Linear Model
		model_lm <- lm(formula=sel_formula, data=lm_list)
		model_f <- summary(model_lm)$fstatistic
		model_c <- summary(model_lm)$coefficients

		# Select variable when there is a significant effect
		if (length(unique(sel_factor)) <= 2 ) {
			if (as.numeric(model_c[2,4]) < pvalue) {
				sel_lm <- c(sel_lm, colnames(feat_matrix)[i])
			}
		} else {
			for (j in 1:length(unique(sel_factor))) {
				if (as.numeric(model_c[j,4]) < pvalue) {
					sel_lm[[as.character(unique(sel_factor)[j])]] <- c(sel_lm[[as.character(unique(sel_factor)[j])]], colnames(feat_matrix)[i])
				}
			}
		}
	}
	
	return(sel_lm)
}



# ---------- Linear Mixed Model with Random Effects ----------
f.select_features_lme <- function(feat_matrix, sel_factor, sel_formula, pvalue=0.005) {
	if (length(unique(sel_factor)) <= 2 ) {
		sel_lme <- NULL
	} else {
		print("Warning. Selection with more than 2 factors is not supported.")
		sel_lme <- list()
	}
	
	for (i in 1:ncol(feat_matrix)) {
		# Create data frame for feature of interest
		lme_list <- data.frame(feature=feat_matrix[,i])
		
		# Linear Model
		#sel_formula=feature ~ age + gender, sel_random_effect=~1|runorder/batch
		model_lme <- lmerTest::lmer(formula=sel_formula, data=lme_list)
		model_t <- summary(model_lme)
		
		# Select variable when there is a significant effect
		if (length(unique(sel_factor)) <= 2 ) {
			if (as.numeric(as.data.frame(model_t$coefficients)[2,'Pr(>|t|)']) < pvalue) {
				sel_lme <- c(sel_lme, colnames(feat_matrix)[i])
			}
		} else {
			for (j in 1:length(unique(sel_factor))) {
				if (as.numeric(as.data.frame(model_t$coefficients)[j,'Pr(>|t|)']) < pvalue) {
					sel_lme[[as.character(unique(sel_factor)[j])]] <- c(sel_lme[[as.character(unique(sel_factor)[j])]], colnames(feat_matrix)[i])
				}
			}
		}
	}
	
	return(sel_lme)
}



# ---------- Feature Selection: Anova + post-hoc Tukey HSD ----------
f.select_features_aov <- function(feat_class, feat_list, conf_level=0.95, pvalue=0.05) {
	sel_aov <- list()
	
	for (i in 1:ncol(feat_list)) {
		# Factor
		aov_fac <- feat_class
		
		# Create data frame for feature of interest
		aov_list <- data.frame(feature=feat_list[,i])
		
		# Attach presence-absence matrix of factors
		aov_list <- cbind(aov_list, as.data.frame(model.matrix(~ 0 + aov_fac)))
		colnames(aov_list)[2:ncol(aov_list)] <- as.character(unique(aov_fac))
		
		# Diagnostic plots
		#hist(aov_list$feature)
		#boxplot(feature ~ aov_fac, data=aov_list, main=as.character(i))
		
		# ANOVA + post-hoc TukeyHSD for each factor
		if (sum(aov_list$feature) > 0) {
			model_aov <- aov(formula=feature ~ aov_fac, data=aov_list)
			model_tukey <- TukeyHSD(x=model_aov, conf.level=conf_level)
			for (j in unique(aov_fac)) {
				if (all(model_tukey[[1]][grepl(j, rownames(model_tukey[[1]])), "p adj"] <= pvalue)) {
					sel_aov[[j]] <- c(sel_aov[[j]], colnames(feat_list)[i])
				}
			}
		}
	}
	
	return(sel_aov)
}



# ---------- LASSO ----------
f.select_features_lasso <- function(feat_matrix, sel_factor, family="gaussian") {
	if (length(unique(sel_factor)) > 2 ) {
		print("Warning. Selection with Lasso only implemented for factor with <= 2 levels.")
	}
	
	# Lasso when alpha=1, Elastic Net when tuning: 0 < alpha < 1
	cv_lasso <- cv.glmnet(x=feat_matrix, y=sel_factor, family=family, alpha=1, nfold=10, type.measure="deviance")
	model_lasso <- glmnet(x=feat_matrix, y=sel_factor, family=family, alpha=1, lambda=cv_lasso$lambda.1se)
	imp_lasso <- coef(model_lasso, s=cv_lasso$lambda.1se)
	sel_lasso <- names(imp_lasso[which(imp_lasso[,1] > 0), ])
	
	return(sel_lasso)
}



# ---------- Elastic Net ----------
f.select_features_elastic_net <- function(feat_matrix, sel_factor, sel_colors, confidence=0.1, tune_length=100, plot_roc_filename=NULL) {
	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}
	
	# Tune alpha and lambda with Lasso and Elastic Net
	model_glmnet <- caret::train(x=as.matrix(feat_matrix), y=sel_factor, method="glmnet",
								 #preProcess=c("center", "scale"),
								 tuneLength=tune_length, trControl=caret::trainControl(method="repeatedcv", number=10, repeats=5, classProbs=TRUE, savePredictions="final"))
	print(paste("Alpha:",as.numeric(model_glmnet$bestTune$alpha)))
	print(paste("Lambda:",as.numeric(model_glmnet$bestTune$lambda)))
	
	# Get variable importances
	imp_glmnet <- varImp(object=model_glmnet)
	rownames(imp_glmnet$importance) <- as.character(rownames(imp_glmnet$importance))
	
	# Names of selected features
	sel_glmnet <- f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_glmnet$importance, confidence=confidence, keepx_min=10)
	
	# Save selected variables
	sel_glmnet[["_selected_variables_"]] <- unique(unlist(f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_glmnet$importance, confidence=confidence, keepx_min=10)))
	
	# Performance measures
	sel_glmnet <- do.call(c, list(sel_glmnet, f.performance_measures_caret(model=model_glmnet, sel_factor=sel_factor, sel_colors=sel_colors)))
	
	dev.off()
	
	return(sel_glmnet)
}



# ---------- PLS ----------
f.select_features_pls <- function(feat_matrix, sel_factor, sel_colors, components=2, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))
	
	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}
	
	# Detach conflicting mixOmics package
	if ("package:mixOmics" %in% search()) detach(package:mixOmics, unload=TRUE)
	
	# Train PLS model
	#tuneGrid=data.frame(ncomp=2)
	model_pls <- caret::train(x=as.matrix(feat_matrix), y=sel_factor, method="pls",
							  preProcess=c("center", "scale"),
							  tuneGrid=data.frame(ncomp=components),
							  tuneLength=tune_length, trControl=caret::trainControl(method="repeatedcv", number=10, repeats=5, classProbs=TRUE, savePredictions="final"))
	print(paste("Number of chosen components:",as.numeric(model_pls$bestTune)))
	
	# Get variable importances
	imp_pls <- varImp(object=model_pls)
	rownames(imp_pls$importance) <- as.character(rownames(imp_pls$importance))
	
	# Names of selected features
	sel_pls <- f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_pls$importance, confidence=quantile_threshold, keepx_min=10)
	
	# Save selected variables
	sel_pls[["_selected_variables_"]] <- unique(unlist(f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_pls$importance, confidence=quantile_threshold, keepx_min=10)))
	
	# Performance measures
	sel_pls <- do.call(c, list(sel_pls, f.performance_measures_caret(model=model_pls, sel_factor=sel_factor, sel_colors=sel_colors)))
	
	dev.off()
	
	return(sel_pls)
}



# ---------- sPLS ----------
f.select_features_spls <- function(feat_matrix, sel_factor, sel_colors, mode="classic", tune_components=2, sel_components=c(1,2), folds_number=10, keepx=c(100,50,25), plot=FALSE, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))
	
	# Set plot off
	if (plot == FALSE) pdf(file=NULL)
	
	# Calculate model
	model_spls <- mixOmics::spls(X=feat_matrix, Y=as.numeric(sel_factor), ncomp=tune_components, mode=mode)
	tune_spls <- mixOmics::perf(object=model_spls, validation="Mfold", folds=folds_number, keepx=keepx, progressBar=TRUE, nrepeat=100)
	
	# Estimate PC axes
	if (plot==TRUE) {
		plot(tune_spls$Q2.total)
		abline(h=0.0975)
		print(tune_spls$Q2.total)
	}
	
	# Plot scores
	if (plot==TRUE) mixOmics::plotIndiv(model_spls, comp=c(1,2), rep.space='X-variate', group=as.numeric(sel_factor), ind.names=sel_factor, legend=TRUE, title="sPLS")
	
	# Partial sPLS
	model_pspls <- mixOmics::spls(X=feat_matrix, Y=as.numeric(sel_factor), ncomp=tune_components, keepX=keepx, mode=mode)
	
	# Plot scores of selected features
	if (plot==TRUE) mixOmics::plotIndiv(model_pspls, comp=c(1,2), rep.space='X-variate', group=as.numeric(sel_factor), ind.names=sel_factor, legend=TRUE, title="partial sPLS")
	
	# Heatmap of selected features
	model_cim <- mixOmics::cim(model_pspls, dist.method=c("euclidean","euclidean"), clust.method=c("complete","ward"), comp=sel_components, xlab="features", ylab="factor", row.names=sel_factor, keysize=c(1,1), keysize.label=0.7, margins=c(4,4), row.cex=0.8, col.cex=0.4)

	# Names of selected features
	if (length(unique(sel_factor)) <= 2 ) {
		sel_spls <- model_cim$col.names
	} else {
		# Generate dendograms
		dend_samp <- as.hclust(model_cim$ddr)
		dend_feat <- as.hclust(model_cim$ddc)
		
		# Generate trees
		tree_samp <- cutree(tree=dend_samp, k=length(unique(sel_factor)))
		tree_feat <- cutree(tree=dend_feat, k=length(unique(sel_factor)))
		
		# Extract features with highest sums in the clusters
		sel_spls <- list()
		mat_sum <- as.data.frame(matrix(nrow=length(unique(tree_samp)), ncol=length(unique(tree_feat))))
		rownames(mat_sum) <- unique(sel_factor)[unique(tree_samp)]
		colnames(mat_sum) <- unique(tree_feat)
		for (i in unique(tree_samp)) {
			for (j in unique(tree_feat)) {
				sel_row <- which(rownames(model_cim$mat) %in% names(tree_samp[tree_samp==i]))
				sel_col <- which(colnames(model_cim$mat) %in% names(tree_feat[tree_feat==j]))
				mat_sum[i,j] <- sum(model_cim$mat[sel_row, sel_col])
			}
			sel_spls[[as.character(unique(sel_factor)[i])]] <- names(tree_feat[tree_feat==which.max(mat_sum[i,])])
		}
		
		sel_spls[["_selected_variables_"]] <- model_cim$col.names
		sel_spls[["_dendrogram_row_"]] <- model_cim$ddr
		sel_spls[["_dendrogram_col_"]] <- model_cim$ddc
	}
	
	# Set plot on
	if (plot == FALSE) dev.off()
	
	# ---------- Performance measures ----------
	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}
	
	# Performance measures, response/observation (CLASS) and prediction (CLASS) vector
	sel_levels <- levels(sel_factor)
	sel_observed <- NULL
	sel_predicted <- NULL
	sel_response <- NULL
	sel_predval <- NULL
	
	# Repeat folds_number times
	for (i in 1:folds_number) {
		# Create training and test datasets for validation
		sample_size <- floor(0.5 * ncol(feat_matrix))
		sample_ind  <- sample(seq_len(ncol(feat_matrix)), size=sample_size)
		
		train_set <- feat_matrix[, sample_ind]
		test_set  <- feat_matrix[, -sample_ind]
		
		# Model
		model_pspls <- mixOmics::spls(X=train_set, Y=as.numeric(sel_factor), ncomp=tune_components, keepX=keepx, mode=mode)
		
		# Prediction
		pred <- predict(model_pspls, newdata=train_set, type="response", decision.values=TRUE, probability=TRUE)

		# Performance measures
		obs <- as.factor(rep(x=sel_factor, times=dim(pred$predict)[3]))
		sel_observed <- c(sel_observed, obs)
		prd <- round(x=as.numeric(pred$predict[,,]), digits=0)
		prd[which(prd<1)] <- 1
		prd[which(prd>length(sel_levels))] <- length(sel_levels)
		sel_predicted <- c(sel_predicted, sel_levels[prd])
		sel_response <- c(sel_response, obs == sel_levels[prd])
		sel_predval <- c(sel_predval, as.numeric(pred$predict[,,]))
	}
	sel_predicted <- as.factor(sel_predicted)
	
	# Plot ROC curves (one vs. rest)
	model_roc <- plot.roc(sel_response, sel_predval, main="Confidence intervals", percent=TRUE, ci=TRUE, print.auc=TRUE)
	model_roc_confidence <- ci.se(model_roc, specificities=seq(0, 100, 1), of="thresholds", thresholds="local maximas") # Confidence
	model_roc_sensitivity <- ci.se(model_roc, specificities=seq(0, 100, 10), of="thresholds", thresholds="local maximas") # Sensitivity
	model_roc_specificity <- ci.sp(model_roc, sensitivities=seq(0, 100, 10), of="thresholds", thresholds="local maximas") # Specificity
	plot(model_roc_confidence, type="shape", col="#1c61b699")
	plot(model_roc_sensitivity, type="bars")
	plot(model_roc_specificity, type="bars")
	
	# Multi-class ROC, response/observation (CLASS) and prediction (0.123) vector
	pred <- pROC::multiclass.roc(response=sel_levels[sel_observed], predictor=sel_predval, levels=unique(sel_factor), percent=FALSE, print.auc=TRUE)
	
	# Multi-class AUC as defined by Hand & Till (2001)
	#print(paste0("Multi-class-AUC: ", round(as.numeric(pred$auc),3)))
	
	# Multi-class Classification rate (= 1 - error rate)
	mcr <- f.classification_rate(sel_factor=sel_observed, predicted=as.numeric(sel_predicted))
	
	# Balanced Error Rate (BER) (= 1 - multi-class classification rate)
	ber <- mlr::measureBER(truth=sel_observed, response=as.numeric(sel_predicted))
	
	# Accuracy
	accuracy <- f.accuracy(sel_factor=sel_observed, predicted=as.numeric(sel_predicted))
	weighted_accuracy <- 0#f.weighted_accuracy(sel_factor=sel_observed, predicted=as.numeric(sel_predicted))
	
	# R2
	#R2 <- f.r2(actual=as.numeric(sel_factor), predicted=as.numeric(pred$predictor))
	R2 <- f.r2(actual=as.numeric(sel_observed), predicted=as.numeric(sel_predicted))
	
	# Precision and Recall Curves
	model_pr <- pr.curve(scores.class0=as.numeric(sel_observed), scores.class1=as.numeric(sel_predicted), curve=TRUE)
	aucpr <- model_pr$auc.integral
	plot(x=model_pr$curve[,1], y=model_pr$curve[,2], xlab="Recall",ylab="Precision", t="l", lwd=2, main="Precision Recall Curve")
	
	# Save metrics
	sel_spls[["_model_roc_"]] <- model_roc
	sel_spls[["_model_roc_confidence_"]] <- model_roc_confidence
	sel_spls[["_model_roc_sensitivity_"]] <- model_roc_sensitivity
	sel_spls[["_model_roc_specificity_"]] <- model_roc_specificity
	sel_spls[["_model_multiclass_roc_"]] <- pred
	sel_spls[["_multiclass_auc_"]] <- as.numeric(pred$auc)
	sel_spls[["_multi_class_rate_"]] <- mcr
	sel_spls[["_balanced_error_rate_"]] <- ber
	sel_spls[["_accuracy_"]] <- accuracy
	sel_spls[["_weighted_accuracy_"]] <- weighted_accuracy
	sel_spls[["_model_r2_"]] <- R2
	sel_spls[["_model_pr_"]] <- model_pr
	sel_spls[["_AUCPR_"]] <- aucpr
	
	# ROC curves for each factor
	pred_fac <- list()
	for (i in sel_levels) {
		#pred_fac[[as.character(i)]] <- as.matrix(model$pred[, c("obs","pred",as.character(i))])
		mat <- data.frame(as.character(sel_levels[sel_observed]), as.character(sel_predicted), as.numeric(sel_predval))
		colnames(mat) <- c("obs", "pred", as.character(i))
		pred_fac[[as.character(i)]] <- mat
		for (j in 1:nrow(pred_fac[[as.character(i)]])) {
			if (pred_fac[[as.character(i)]][j, "obs"] == as.character(i)) {
				if (pred_fac[[as.character(i)]][j, "pred"] == as.character(i)) {
					pred_fac[[as.character(i)]][j, "pred"] <- 1
				} else {
					pred_fac[[as.character(i)]][j, "pred"] <- 0
				}
				pred_fac[[as.character(i)]][j, "obs"] <- 1
			} else {
				if (pred_fac[[as.character(i)]][j, "pred"] != as.character(i)) {
					pred_fac[[as.character(i)]][j, "pred"] <- 1
				} else {
					pred_fac[[as.character(i)]][j, "pred"] <- 0
				}
				pred_fac[[as.character(i)]][j, "obs"] <- 0
			}
		}
	}
	
	model_roc <- list()
	model_roc_confidence <- list()
	model_roc_sensitivity <- list()
	model_roc_specificity <- list()
	for (i in 1:length(sel_levels)) {
		if (i==1) add=FALSE else add=TRUE
		model_roc[[as.character(sel_levels[i])]] <- plot.roc(as.numeric(pred_fac[[sel_levels[i]]][,"obs"]), as.numeric(pred_fac[[sel_levels[i]]][,as.character(sel_levels[i])]), main="ROC curves of levels", col=sel_colors[i], percent=TRUE, ci=TRUE, print.auc=FALSE, add=add)
		model_roc_confidence[[as.character(sel_levels[i])]] <- ci.se(model_roc[[as.character(sel_levels[i])]], specificities=seq(0, 100, 1), of="thresholds", thresholds="local maximas") # Confidence
		#model_roc_sensitivity[[as.character(sel_levels[i])]] <- ci.se(model_roc[[as.character(sel_levels[i])]], specificities=seq(0, 100, 10), of="thresholds", thresholds="local maximas") # Sensitivity
		#model_roc_specificity[[as.character(sel_levels[i])]] <- ci.sp(model_roc[[as.character(sel_levels[i])]], sensitivities=seq(0, 100, 10), of="thresholds", thresholds="local maximas") # Specificity
		plot(model_roc_confidence[[as.character(sel_levels[i])]], type="shape", col=rgb(col2rgb(sel_colors[i])[1]/255, col2rgb(sel_colors[i])[2]/255, col2rgb(sel_colors[i])[3]/255, alpha=0.3), lty=0, no.roc=TRUE)
		#plot(model_roc_sensitivity[[as.character(sel_levels[i])]], type="bars")
		#plot(model_roc_specificity[[as.character(sel_levels[i])]], type="bars")
	}
	legend("bottomright", legend=paste0(sel_levels, ", AUC: ", round(as.numeric(unlist(lapply(model_roc, function(x) x$auc))), 1), "%"), col=sel_colors, lwd=2, cex=0.75)
	
	# Precision and Recall Curves
	model_pr <- list()
	plot(xlim=c(0,1), ylim=c(0,1), x=NULL, y=NULL, xlab="Recall", ylab="Precision", main="Precision Recall Curves of levels")
	for (i in 1:length(sel_levels)) {
		model_pr[[as.character(sel_levels[i])]] <- pr.curve(scores.class0=as.numeric(pred_fac[[sel_levels[i]]][,"obs"]), scores.class1=as.numeric(pred_fac[[sel_levels[i]]][,as.character(sel_levels[i])]), curve=TRUE)
		lines(x=model_pr[[as.character(sel_levels[i])]]$curve[,1], y=model_pr[[as.character(sel_levels[i])]]$curve[,2], xlab="Recall",ylab="Precision", t="l", col=sel_colors[i], lwd=2)
	}
	legend("topright", legend=paste0(sel_levels, ", AUC-PR: ", round(as.numeric(unlist(lapply(model_pr, function(x) x$auc.integral*100))), 1), "%"), col=sel_colors, lwd=2, cex=0.75)
	
	# Save metrics
	sel_spls[["_model_factor_roc_"]] <- model_roc
	sel_spls[["_model_factor_roc_confidence_"]] <- model_roc_confidence
	sel_spls[["_model_factor_roc_sensitivity_"]] <- model_roc_sensitivity
	sel_spls[["_model_factor_roc_specificity_"]] <- model_roc_specificity
	sel_spls[["_model_factor_pr_"]] <- model_pr
	
	dev.off()
	
	return(sel_spls)
}



# ---------- sPLS-DA ----------
f.select_features_splsda <- function(feat_matrix, sel_factor, sel_colors, tune_components=2, sel_components=c(1,2), folds_number=10, keepx=c(100,50,25), plot=FALSE, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))
	
	# Set plot off
	if (plot == FALSE) pdf(file=NULL)
	
	# Calculate model
	model_splsda <- mixOmics::splsda(X=feat_matrix, Y=sel_factor, ncomp=tune_components)
	tune_splsda <- mixOmics::perf(object=model_splsda, validation="Mfold", folds=folds_number, progressBar=TRUE, nrepeat=10, cpus=1)
	
	# Estimate PC axes
	#if (plot==TRUE) {
	#	plot(tune_splsda$Q2.total)
	#	abline(h=0.0975)
	#	print(tune_splsda$Q2.total)
	#}
	
	# Plot scores
	if (plot==TRUE) mixOmics::plotIndiv(model_splsda, comp=c(1,2), rep.space='X-variate', group=as.numeric(sel_factor), ind.names=sel_factor, legend=TRUE, title="sPLS-DA")
	
	# Partial sPLS-DA
	model_psplsda <- mixOmics::splsda(X=feat_matrix, Y=as.numeric(sel_factor), ncomp=tune_components, keepX=keepx)
	tune_psplsda <- mixOmics::perf(object=model_psplsda, validation="Mfold", folds=folds_number, progressBar=TRUE, nrepeat=10, cpus=1)
	
	# Plot scores of selected features
	if (plot==TRUE) mixOmics::plotIndiv(model_psplsda, comp=c(1,2), rep.space='X-variate', group=as.numeric(sel_factor), ind.names=sel_factor, legend=TRUE, title="partial sPLS-DA")
	
	# Heatmap of selected features
	model_cim <- mixOmics::cim(model_psplsda, dist.method=c("euclidean","euclidean"), clust.method=c("complete","ward"), comp=sel_components, xlab="features", ylab="factor", row.names=sel_factor, keysize=c(1,1), keysize.label=0.7, margins=c(4,4), row.cex=0.8, col.cex=0.4)
	
	# Generate dendograms
	dend_samp <- as.hclust(model_cim$ddr)
	dend_feat <- as.hclust(model_cim$ddc)
	
	# Generate trees
	tree_samp <- cutree(tree=dend_samp, k=length(unique(sel_factor)))
	tree_feat <- cutree(tree=dend_feat, k=length(unique(sel_factor)))
	
	# Extract features with highest sums in the clusters
	sel_splsda <- list()
	mat_sum <- as.data.frame(matrix(nrow=length(unique(tree_samp)), ncol=length(unique(tree_feat))))
	rownames(mat_sum) <- unique(sel_factor)[unique(tree_samp)]
	colnames(mat_sum) <- unique(tree_feat)
	for (i in unique(tree_samp)) {
		for (j in unique(tree_feat)) {
			sel_row <- which(rownames(model_cim$mat) %in% names(tree_samp[tree_samp==i]))
			sel_col <- which(colnames(model_cim$mat) %in% names(tree_feat[tree_feat==j]))
			mat_sum[i,j] <- sum(model_cim$mat[sel_row, sel_col])
		}
		sel_splsda[[as.character(unique(sel_factor)[i])]] <- names(tree_feat[tree_feat==which.max(mat_sum[i,])])
	}
	
	sel_splsda[["_selected_variables_"]] <- model_cim$col.names
	sel_splsda[["_dendrogram_row_"]] <- model_cim$ddr
	sel_splsda[["_dendrogram_col_"]] <- model_cim$ddc

	# Set plot on
	if (plot == FALSE) dev.off()
	
	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}
	
	# Seletion levels
	sel_levels <- levels(as.factor(sel_factor))
	
	# Create prediction matrix
	model_pred <- NULL
	for (i in 1:dim(tune_psplsda$predict[[length(tune_psplsda$predict)]])[3]) model_pred <- rbind(model_pred, tune_psplsda$predict[[length(tune_psplsda$predict)]][,,i])
	colnames(model_pred) <- sel_levels
	rownames(model_pred) <- rep(as.character(sel_factor), times=dim(tune_psplsda$predict[[length(tune_psplsda$predict)]])[3])
	model_pred <- scales::rescale(x=model_pred, to=c(0,1))
	colnames(model_pred) <- paste0(colnames(model_pred), "_pred_")
	pred_true <- data.frame(dummies::dummy(gsub(x=rownames(model_pred), pattern='\\..*', replacement='')))
	colnames(pred_true) <- paste0(sel_levels, "_true")
	model_pred <- cbind(pred_true, model_pred)
	
	# Create classification objects
	#sel_obs <- NULL; for (i in which(((model_pred[,c((nlevels(sel_levels)+1):ncol(model_pred))] * (model_pred[,c(1:(nlevels(sel_levels)))]==1)) > 0))) sel_obs <- c(sel_obs, gsub(x=names(unlist(model_pred[,c(1:(nlevels(sel_levels)))]))[i], pattern='_.*', replacement=''))
	sel_obs <- gsub(x=rownames(model_pred), pattern='\\.\\d.*', replacement='')
	sel_prob <- model_pred[,c((nlevels(sel_levels)+1):ncol(model_pred))][((model_pred[,c((nlevels(sel_levels)+1):ncol(model_pred))] * (model_pred[,c(1:(nlevels(sel_levels)))]==1)) > 0)]
	#sel_pred <- NULL; for (i in c(1:dim(tune_psplsda$class$max.dist)[2])) { sel_pred <- c(sel_pred, as.character(sel_levels[as.numeric(apply(X=tune_psplsda$class$max.dist[,i,], MARGIN=1, FUN=function(x) { names(which.max(table(x))) }))])) }
	#sel_pred <- NULL; for (i in c(1:dim(tune_splsda$class$max.dist)[2])) { sel_pred <- c(sel_pred, as.character(tune_splsda$class$max.dist[,i,dim(tune_splsda$class$max.dist)[3]])) }
	sel_pred <- as.character(apply(X=model_pred[,c((length(sel_levels)+1):ncol(model_pred))], MARGIN=1, FUN=function(x) { gsub(x=names(which.max(x)), pattern='_.*', replacement='') }))
	
	# Build caret like object
	model_caret <- list()
	model_caret$pred <- data.frame(pred=sel_pred, obs=sel_obs)
	model_caret$pred <- cbind(model_caret$pred, model_pred[,c((length(sel_levels)+1):ncol(model_pred))])
	colnames(model_caret$pred) <- c("pred", "obs", as.character(sel_levels))
	
	# Performance measures
	sel_splsda <- do.call(c, list(sel_splsda, f.performance_measures_caret(model=model_caret, sel_factor=sel_factor, sel_colors=sel_colors)))
	
	dev.off()

	return(sel_splsda)
}



# ---------- PLS-DA ----------
f.select_features_plsda <- function(feat_matrix, sel_factor, sel_colors, tune_components=2, sel_components=c(1,2), folds_number=10, dist.method=c("euclidean","euclidean"), clust.method=c("complete","ward"), plot=FALSE, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))
	
	# Set plot off
	if (plot == FALSE) pdf(file=NULL)
	
	# Calculate model
	model_plsda <- mixOmics::plsda(X=feat_matrix, Y=sel_factor, ncomp=tune_components)
	tune_plsda <- mixOmics::perf(object=model_plsda, validation="Mfold", folds=folds_number, progressBar=TRUE, nrepeat=10, cpus=1)
	
	# Estimate PC axes
	#if (plot==TRUE) {
	#	plot(tune_plsda$Q2.total)
	#	abline(h=0.0975)
	#	print(tune_plsda$Q2.total)
	#}
	
	# Plot scores
	if (plot==TRUE) mixOmics::plotIndiv(model_plsda, comp=c(1,2), rep.space='X-variate', group=as.numeric(sel_factor), ind.names=sel_factor, legend=TRUE, title="pls-DA")
	
	# Partial PLS-DA
	model_pplsda <- mixOmics::plsda(X=feat_matrix, Y=as.numeric(sel_factor), ncomp=tune_components)
	tune_pplsda <- mixOmics::perf(object=model_pplsda, validation="Mfold", folds=folds_number, progressBar=TRUE, nrepeat=10, cpus=1)
	
	# Plot scores of selected features
	if (plot==TRUE) mixOmics::plotIndiv(model_pplsda, comp=c(1,2), rep.space='X-variate', group=as.numeric(sel_factor), ind.names=sel_factor, legend=TRUE, title="partial pls-DA")
	
	# Heatmap of selected features
	model_cim <- mixOmics::cim(model_pplsda, dist.method=dist.method, clust.method=clust.method, comp=sel_components, xlab="features", ylab="factor", row.names=sel_factor, keysize=c(1,1), keysize.label=0.7, margins=c(4,4), row.cex=0.8, col.cex=0.4)
	
	# Generate dendograms
	dend_samp <- as.hclust(model_cim$ddr)
	dend_feat <- as.hclust(model_cim$ddc)
	
	# Generate trees
	tree_samp <- cutree(tree=dend_samp, k=length(unique(sel_factor)))
	tree_feat <- cutree(tree=dend_feat, k=length(unique(sel_factor)))
	
	# Extract features with highest sums in the clusters
	sel_plsda <- list()
	mat_sum <- as.data.frame(matrix(nrow=length(unique(tree_samp)), ncol=length(unique(tree_feat))))
	rownames(mat_sum) <- unique(sel_factor)[unique(tree_samp)]
	colnames(mat_sum) <- unique(tree_feat)
	for (i in unique(tree_samp)) {
		for (j in unique(tree_feat)) {
			sel_row <- which(rownames(model_cim$mat) %in% names(tree_samp[tree_samp==i]))
			sel_col <- which(colnames(model_cim$mat) %in% names(tree_feat[tree_feat==j]))
			mat_sum[i,j] <- sum(model_cim$mat[sel_row, sel_col])
		}
		sel_plsda[[as.character(unique(sel_factor)[i])]] <- names(tree_feat[tree_feat==which.max(mat_sum[i,])])
	}
	
	sel_plsda[["_selected_variables_"]] <- model_cim$col.names
	sel_plsda[["_dendrogram_row_"]] <- model_cim$ddr
	sel_plsda[["_dendrogram_col_"]] <- model_cim$ddc

	# Set plot on
	if (plot == FALSE) dev.off()
	
	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}
	
	# Seletion levels
	sel_levels <- levels(as.factor(sel_factor))
	
	# Create prediction matrix
	model_pred <- NULL
	for (i in 1:dim(tune_pplsda$predict[[length(tune_pplsda$predict)]])[3]) model_pred <- rbind(model_pred, tune_pplsda$predict[[length(tune_pplsda$predict)]][,,i])
	colnames(model_pred) <- sel_levels
	rownames(model_pred) <- rep(as.character(sel_factor), times=dim(tune_pplsda$predict[[length(tune_pplsda$predict)]])[3])
	model_pred <- scales::rescale(x=model_pred, to=c(0,1))
	colnames(model_pred) <- paste0(colnames(model_pred), "_pred_")
	pred_true <- data.frame(dummies::dummy(gsub(x=rownames(model_pred), pattern='\\..*', replacement='')))
	colnames(pred_true) <- paste0(sel_levels, "_true")
	model_pred <- cbind(pred_true, model_pred)
	
	# Create classification objects
	#sel_obs <- NULL; for (i in which(((model_pred[,c((nlevels(sel_levels)+1):ncol(model_pred))] * (model_pred[,c(1:(nlevels(sel_levels)))]==1)) > 0))) sel_obs <- c(sel_obs, gsub(x=names(unlist(model_pred[,c(1:(nlevels(sel_levels)))]))[i], pattern='_.*', replacement=''))
	sel_obs <- gsub(x=rownames(model_pred), pattern='\\..*', replacement='')
	sel_prob <- model_pred[,c((nlevels(sel_levels)+1):ncol(model_pred))][((model_pred[,c((nlevels(sel_levels)+1):ncol(model_pred))] * (model_pred[,c(1:(nlevels(sel_levels)))]==1)) > 0)]
	#sel_pred <- NULL; for (i in c(1:dim(tune_pplsda$class$max.dist)[2])) { sel_pred <- c(sel_pred, as.character(sel_levels[as.numeric(apply(X=tune_pplsda$class$max.dist[,i,], MARGIN=1, FUN=function(x) { names(which.max(table(x))) }))])) }
	#sel_pred <- NULL; for (i in c(1:dim(tune_plsda$class$max.dist)[2])) { sel_pred <- c(sel_pred, as.character(tune_plsda$class$max.dist[,i,dim(tune_plsda$class$max.dist)[3]])) }
	sel_pred <- as.character(apply(X=model_pred[,c((length(sel_levels)+1):ncol(model_pred))], MARGIN=1, FUN=function(x) { gsub(x=names(which.max(x)), pattern='_.*', replacement='') }))
	
	# Build caret like object
	model_caret <- list()
	model_caret$pred <- data.frame(pred=sel_pred, obs=sel_obs)
	model_caret$pred <- cbind(model_caret$pred, model_pred[,c((length(sel_levels)+1):ncol(model_pred))])
	colnames(model_caret$pred) <- c("pred", "obs", as.character(sel_levels))
	
	# Performance measures
##	sel_plsda <- do.call(c, list(sel_plsda, f.performance_measures_caret(model=model_caret, sel_factor=sel_factor, sel_colors=sel_colors)))
	
	dev.off()

	return(sel_plsda)
}



# ---------- PLS-DA CARET ----------
f.select_features_plsda_caret <- function(feat_matrix, sel_factor, sel_colors, ncomp=10, quantile_threshold, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))
	
	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}
	
	# Train PLS-DA model
	model_plsda <- caret::train(x=as.matrix(feat_matrix), y=sel_factor, method="plsda", ncomp=ncomp, probMethod="Bayes",
							    trControl=caret::trainControl(method="repeatedcv", number=10, repeats=5, classProbs=TRUE, savePredictions="final"))
	
	# Get variable importances
	imp_rf <- varImp(object=model_rf)
	rownames(imp_rf$importance) <- as.character(rownames(imp_rf$importance))
	
	# Save names of selected features
	sel_rf <- f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_rf$importance, confidence=quantile_threshold, keepx_min=10)
	
	# Save selected variables
	sel_rf[["_selected_variables_"]] <- unique(unlist(f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_rf$importance, confidence=quantile_threshold, keepx_min=10)))
	
	# Performance measures
	sel_rf <- do.call(c, list(sel_rf, f.performance_measures_caret(model=model_rf, sel_factor=sel_factor, sel_colors=sel_colors)))
	
	dev.off()
	
	return(sel_rf)
}



# ---------- oPLS ----------
f.select_features_opls <- function(feat_matrix, sel_factor, chosen_features) {
	# Model
	model_opls <- ropls::opls(x=feat_matrix, y=sel_factor, scaleC="standard", predI=1, orthoI=1)
	ropls::plot(x=model_opls, y=sel_factor, typeVc="x-score", parLabVc=as.character(sel_factor), parDevNewL=FALSE, parLayL=FALSE, parCexN=0.5, parEllipsesL=FALSE, parTitleL=TRUE)
	sel_opls <- names(head(model_opls@vipVn[order(model_opls@vipVn, decreasing=TRUE)], n=chosen_features))
	
	return(sel_opls)
}



# ---------- Random Forest ----------
f.select_features_random_forest <- function(feat_matrix, sel_factor, sel_colors, tune_length, quantile_threshold, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))
	
	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}
	
	# Train RF model
	model_rf <- caret::train(x=as.matrix(feat_matrix), y=sel_factor, method="rf", importance=TRUE, proximity=TRUE,
							 tuneLength=tune_length, trControl=caret::trainControl(method="repeatedcv", number=10, repeats=5, classProbs=TRUE, savePredictions="final"))
	
	# Get variable importances
	imp_rf <- varImp(object=model_rf)
	rownames(imp_rf$importance) <- as.character(rownames(imp_rf$importance))
	
	# Save names of selected features
	sel_rf <- f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_rf$importance, confidence=quantile_threshold, keepx_min=10)
	
	# Save selected variables
	sel_rf[["_selected_variables_"]] <- unique(unlist(f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_rf$importance, confidence=quantile_threshold, keepx_min=10)))
	
	# Performance measures
	sel_rf <- do.call(c, list(sel_rf, f.performance_measures_caret(model=model_rf, sel_factor=sel_factor, sel_colors=sel_colors)))
	
	dev.off()
	
	return(sel_rf)
}



# ---------- Random Forest using Boruta selection method ----------
f.select_features_random_forest_boruta <- function(feat_matrix, sel_factor, sel_colors, sel_method="best", tune_length=10, pvalue=0.05, mcAdj=TRUE, maxRuns=1000, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))
	
	#library(Boruta)
	#library(rpart)
	
	#sel_boruta <- list()
	#for (i in 1:tune_length) {
	sel_boruta <- foreach(i=c(1:tune_length)) %dopar% {
		sel_boruta_rf <- list()
		# Perform Random Forest with default p-value and Bonferroni multiple comparisons adjustment
		model_sel_boruta <- Boruta::Boruta(x=as.matrix(feat_matrix), y=sel_factor, pValue=pvalue, mcAdj=mcAdj, maxRuns=maxRuns, doTrace=0, holdHistory=TRUE, getImp=Boruta::getImpRfZ)
		
		#getSelectedAttributes(model_sel_boruta, withTentative=FALSE)
		#model_sel_boruta <- attStats(model_sel_boruta)
		
		# Names of selected features
		sel_boruta_rf[["model"]] <- model_sel_boruta
		sel_boruta_rf[["_selected_variables_"]] <- names(model_sel_boruta$finalDecision[which(model_sel_boruta$finalDecision=="Confirmed")])
		
		# Build regression tree
		model_rtree <- rpart::rpart(formula=as.formula(paste0("sel_factor ~ ", paste(sel_boruta_rf[["_selected_variables_"]], collapse=' + '))), data=as.data.frame(feat_matrix), model=TRUE)
		sel_boruta_rf[["_model_rtree_"]] <- model_rtree
		
		# Calculate actual vs. predicted
		actual <- sel_factor
		sel_boruta_rf[["_actual_"]] <- actual
		actual <- as.numeric(actual)
		
		sel_boruta_rf[["_class_probabilities_"]] <- as.data.frame(rpart:::predict.rpart(object=model_rtree, type="prob"))
		sel_boruta_rf[["_predicted_"]] <- rpart:::predict.rpart(object=model_rtree, type="class")
		predicted <- rpart:::predict.rpart(object=model_rtree, type="vector")
		predicted <- as.numeric(predicted)
		
		# R-squared
		print(paste(i, caret::postResample(actual, predicted)))
		sel_boruta_rf[["_R2_"]] <- caret::postResample(actual, predicted)[2]
		
		# Root Mean Square Error
		sel_boruta_rf[["_RMSE_"]] <- caret::postResample(actual, predicted)[1]
		
		# Mean Absolute Error
		sel_boruta_rf[["_MAE_"]] <- caret::postResample(actual, predicted)[3]
		
		#sel_boruta <- c(sel_boruta, sel_boruta_rf)
		return(sel_boruta_rf)
	}

	# Save metrics
	#sel_boruta[["_R2_"]] <- as.numeric(subListExtract(sel_boruta, '_R2_', simplify=TRUE))
	#sel_boruta[["_RMSE_"]] <- as.numeric(subListExtract(sel_boruta, '_RMSE_', simplify=TRUE))
	#sel_boruta[["_MAE_"]] <- as.numeric(subListExtract(sel_boruta, '_MAE_', simplify=TRUE))
	
	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}

	# Select best performing model
	if (sel_method == "best") {
		# Best model
		sel_vars <- as.numeric(subListExtract(sel_boruta, '_R2_', simplify=TRUE))
		best_model <- which(sel_vars==max(sel_vars))[1]
		sel_vars <- sel_boruta[[best_model]][['_selected_variables_']]
		
		# Selection levels
		sel_levels <- levels(as.factor(sel_factor))
		
		# Build data frame
		#model_pred <- data.frame(mtry=rep(1:tune_length, each=length(sel_factor)))
		#model_pred$pred <- unlist(subListExtract(sel_boruta, '_actual_', simplify=FALSE))
		#model_pred$obs <- unlist(subListExtract(sel_boruta, '_predicted_', simplify=FALSE))
		#model_pred <- cbind(model_pred, do.call(rbind, subListExtract(sel_boruta, '_class_probabilities_', simplify=FALSE)))
		
		model_pred <- data.frame(mtry=rep(best_model, each=length(sel_factor)))
		model_pred$pred <- sel_boruta[[best_model]][['_predicted_']]
		model_pred$obs <- sel_factor
		model_pred <- cbind(model_pred, sel_boruta[[best_model]][['_class_probabilities_']])
		rownames(model_pred) <- make.unique(as.character(sel_factor))
		
		# Create classification objects
		sel_obs <- model_pred$obs
		sel_prob <- model_pred[,c((nlevels(sel_levels)+1):ncol(model_pred))][((model_pred[,c((nlevels(sel_levels)+1):ncol(model_pred))] * (model_pred[,c(1:(nlevels(sel_levels)))]==1)) > 0)]
		sel_pred <- as.character(apply(X=model_pred[,c((length(sel_levels)+1):ncol(model_pred))], MARGIN=1, FUN=function(x) { gsub(x=names(which.max(x)), pattern='_.*', replacement='') }))
		
		# Build caret like object
		model_caret <- list()
		model_caret$pred <- data.frame(pred=sel_pred, obs=sel_obs)
		model_caret$pred <- cbind(model_caret$pred, model_pred[,c(4:ncol(model_pred))])
		colnames(model_caret$pred) <- c("pred", "obs", as.character(sel_levels))
		
		sel_boruta[['_R2_']] <- sel_boruta[[best_model]][['_R2_']]
		sel_boruta[['_RMSE_']] <- sel_boruta[[best_model]][['_RMSE_']]
		sel_boruta[['_MAE_']] <- sel_boruta[[best_model]][['_MAE_']]
		
		sel_boruta <- do.call(c, list(sel_boruta, f.performance_measures_caret(model=model_caret, sel_factor=sel_factor, sel_colors=sel_colors)))
	# Select common variables
	} else if (sel_method == "common") {
		print("Not implemented.")
		return(NA)
	# Select all variables
	} else {
		# All variables
		sel_vars <- unlist(subListExtract(sel_boruta, '_selected_variables_', simplify=FALSE))
		sel_vars <- sel_vars[! duplicated(sel_vars)]

		# Selection levels
		sel_levels <- levels(as.factor(sel_factor))
		
		# Build regression tree
		model_rtree <- rpart::rpart(formula=as.formula(paste0("sel_factor ~ ", paste(sel_vars, collapse=' + '))), data=as.data.frame(feat_matrix), model=TRUE)

		# Build data frame
		model_pred <- data.frame(mtry=rep(1, each=length(sel_factor)))
		model_pred$pred <- as.numeric(sel_factor)
		model_pred$obs <- rpart:::predict.rpart(object=model_rtree, type="vector")
		model_pred <- cbind(model_pred, as.data.frame(rpart:::predict.rpart(object=model_rtree, type="prob")))
		
		# Create classification objects
		sel_obs <- sel_levels[rpart:::predict.rpart(object=model_rtree, type="vector")]
		sel_prob <- as.data.frame(rpart:::predict.rpart(object=model_rtree, type="prob"))
		sel_pred <- as.character(sel_factor)
		
		# Build caret like object
		model_caret <- list()
		model_caret$pred <- data.frame(pred=sel_pred, obs=sel_obs)
		model_caret$pred <- cbind(model_caret$pred, model_pred[,c(4:ncol(model_pred))])
		colnames(model_caret$pred) <- c("pred", "obs", as.character(sel_levels))
		
		#sel_boruta <- do.call(c, list(sel_boruta, f.performance_measures_caret(model=model_caret, sel_factor=sel_factor, sel_colors=sel_colors)))
	}
	
	sel_boruta[["_selected_variables_"]] <- as.character(sel_vars)
	
	dev.off()
	
	return(sel_boruta)
}



# ---------- Support Vector Machines ----------
f.select_features_svm <- function(feat_matrix, sel_factor, sel_colors, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL) {
	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}
	
	# Train SVM model
	model_svm <- caret::train(x=as.matrix(feat_matrix), y=as.factor(sel_factor), method="svmLinear2", shrinking=TRUE, probability=TRUE, fitted=TRUE,
							  tuneLength=tune_length, trControl=caret::trainControl(method="repeatedcv", number=10, repeats=5, classProbs=TRUE, savePredictions="final"))
	model_svmi <- varImp(object=model_svm) #RFE Recursive Feature Elimination
	
	# Get variable importances
	rownames(model_svmi$importance) <- as.character(rownames(model_svmi$importance))
	sel_svm <- f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=model_svmi$importance, confidence=quantile_threshold, keepx_min=10)
	
	# Save selected variables
	sel_svm[["_selected_variables_"]] <- unique(unlist(f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=model_svmi$importance, confidence=quantile_threshold, keepx_min=10)))
	
	# Performance measures
	sel_svm <- do.call(c, list(sel_svm, f.performance_measures_caret(model=model_svm, sel_factor=sel_factor, sel_colors=sel_colors)))
	
	dev.off()
	
	return(sel_svm)
}



# ---------- Self Organizing Maps ----------
f.select_features_som <- function(feat_matrix, sel_factor, sel_colors, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL) {
	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}
	
	# Train SOM model
	model_som <- caret::train(x=as.matrix(feat_matrix), y=as.factor(sel_factor), method="xyf",
							  tuneLength=tune_length, trControl=caret::trainControl(method="repeatedcv", number=10, repeats=5, classProbs=TRUE, savePredictions="final"))
	model_somi <- varImp(object=model_som)
	
	# Get variable importances
	rownames(model_somi$importance) <- as.character(rownames(model_somi$importance))
	sel_som <- f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=model_somi$importance, confidence=quantile_threshold, keepx_min=10)
	
	# Save selected variables
	sel_som[["_selected_variables_"]] <- unique(unlist(f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=model_somi$importance, confidence=quantile_threshold, keepx_min=10)))
	
	# Performance measures
	sel_som <- do.call(c, list(sel_som, f.performance_measures_caret(model=model_som, sel_factor=sel_factor, sel_colors=sel_colors)))
	
	dev.off()
	
	return(sel_som)
}



# ---------- Select features from model ----------
f.select_features_from_model <- function(feat_list, feat_class, model_varimp, keepx_min, confidence=0.95) {
	if (ncol(as.data.frame(model_varimp)) < 2 ) {
		model_varimp <- as.data.frame(model_varimp)
		model_varimp <- cbind(model_varimp, model_varimp[,1])
		colnames(model_varimp) <- as.character(unique(feat_class))
	}
	
	# Create list
	sel_list <- list()
	
	# Limit selected features between keepx_max and keepx_min
	for (i in unique(feat_class)) {
		elements <- colnames(feat_list)[which(model_varimp[,i] >= confidence * max(model_varimp[,i]))]
		
		#if (length(elements) > keepx_max) {
		#  sel_list[[i]] <- colnames(feat_list)[base::order(model_varimp[, i], decreasing=TRUE)[1:keepx_max]]
		#} else
		if (length(elements) < keepx_min) {
			sel_list[[i]] <- colnames(feat_list)[base::order(model_varimp[, i], decreasing=TRUE)[1:keepx_min]]
		} else {
			sel_list[[i]] <- elements
		}
		
		sel_list[[i]] <- sort(sel_list[[i]])
	}
	
	# Return selected features
	return(sel_list)
}



# ---------- Draw heatmap of selected features ----------
f.count.selected_features <- function(sel_feat) {
	# Return number of selected features
	return(length(unique(unlist(sel_feat))))
}



# ---------- Draw heatmap of selected features ----------
f.heatmap.selected_features <- function(feat_list, sel_feat, sel_names=NULL, sample_colors=NULL, filename, main, scale="col", plot_width=6, plot_height=5, cex_col=0.5, cex_row=0.7) {
	# Use existing dendrogram for rows
	if (any(names(sel_feat) %in% '_dendrogram_row_')) {
		print("Using existing dendrogram for clustering of rows.")
		rowv = sel_feat$'_dendrogram_row_'
		sel_feat$'_dendrogram_row_' <- NULL
	} else {
		rowv = NULL
	}
	
	# Use existing dendrogram for columns
	if (any(names(sel_feat) %in% '_dendrogram_col_')) {
		print("Using existing dendrogram for clustering of columns.")
		colv = sel_feat$'_dendrogram_col_'
		sel_feat$'_dendrogram_col_' <- NULL
	} else {
		colv = NULL
	}
	
	# Plot '_selected_variables_' if they exist
	if (any(names(sel_feat) %in% '_selected_variables_')) {
		print("Using matrix with selected variables.")
		sel_list <- scale((feat_list[, which(colnames(feat_list) %in% sel_feat[["_selected_variables_"]])]),scale=T,center=T)
	} else {
		# Data frame with only selected features
		sel_list <- as.data.frame(feat_list[,c(sort(as.character(unique(unlist(sel_feat))))), drop=FALSE])
	}
	
	# Use sel_names
	if (! is.null(sel_names)) {
		colnames(sel_list) <- sel_names
	}
	
	# Use colors for samples (row)
	if (is.null(sample_colors)) {
		sample_colors <- "black"
	}
	
	# Clustering of rows and columns
	if (is.null(rowv)) {
		rowv = as.dendrogram(hclust(dist(scale(sel_list),method="euclidean"),method="complete"), center=T)
		#rowv = as.dendrogram(hclust(dist(as.numeric(as.factor(mymetadata$SpecCode))),method="ward.D"), center=T), offsetRow=0, #<< samples need to be sorted by SpecCode (not clustered)
	}
	if (is.null(colv)) {
		colv = as.dendrogram(hclust(dist(t(scale(sel_list)),method="euclidean"),method="ward.D"), center=T)
	}
	
	# Draw heatmap
	if (! is.null(filename)) pdf(file=as.character(filename), encoding="ISOLatin1", pointsize=10, width=plot_width, height=plot_height, family="Helvetica")
	heatmap.2(x=as.matrix(sel_list), scale=scale, cexRow=cex_row, cexCol=cex_col, main=main,
			  Rowv=rowv, offsetRow=0, colRow=sample_colors,
			  Colv=colv, offsetCol=0,
			  col=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256),
			  #trace="none", margins=c(2,4),
			  trace="none", margins=c(max(nchar(colnames(sel_list)))/5, max(nchar(rownames(sel_list)))/3),
			  key=TRUE, key.title="Color key", density.info='density', denscol="black")
	if (! is.null(filename)) dev.off()
}



# ---------- Perform Dynamic Time Warping ----------
f.dynamic_time_warping <- function(loess_list, sel_factor, loess_smoothing=0.9, family="gaussian", clust_type="hierarchical", clust_method="average", curve_types=8, distance="dtw", plot_filename=NULL) {
	#library(dtwclust)
	
	# Loess smoothing for each feature
	for (i in 1:ncol(loess_list)) {
		# Loess smoothing parameter
		model_opt_smooth <- NULL
		model_opt_smooth$par <- loess_smoothing
		
		# Loess curve
		model_loess <- loess(loess_list[,i] ~ as.numeric(as.factor(sel_factor)), span=model_opt_smooth$par, parametric=FALSE, normalize=TRUE, family=family, control=loess.control(surface="direct"))
		
		# Save smoothed values
		loess_list[,i] <- predict(model_loess)
	}
	
	# Dynamic Time Warping Zscores
	dtw_list <- list()
	for (i in 1:ncol(loess_list)) {
		dtw_list[[i]] <- loess_list[,i]
	}
	dtw_list <- zscore(dtw_list)
	
	# Hierarchical clustering of similar curves
	hc_sbd <- tsclust(dtw_list, type=clust_type, k=curve_types,
					  preproc=zscore,
					  distance=distance, centroid=shape_extraction,
					  control=hierarchical_control(method=clust_method))
	
	# Plot
	if (! is.null(plot_filename)) {
		pdf(file=plot_filename, encoding="ISOLatin1", pointsize=10, width=10, height=10, family="Helvetica")
		plot(hc_sbd, type="sc")
		plot(hc_sbd, type="centroids")
		dev.off()
	}
	
	return(hc_sbd)
}



# ---------- Sunburst plot ----------
sunBurstPlotFromSubstanceClasses <- function(classifierClasses, numberOfSpectra, colorStart = 0.2, colorAlpha = 0.5){
	level <- unlist(lapply(X = strsplit(x = classifierClasses, split = "; "), FUN = length))
	
	## all class levels
	classesAndSubClasses <- lapply(X = strsplit(x = classifierClasses, split = "; "), FUN = function(x){
		sapply(X = seq_along(x), FUN = function(y){paste(x[1:y], collapse = "; ")})
	})
	classesByLevel <- list()
	labelsByLevel <- list()
	for(levelHere in seq_len(max(level))){
		classesByLevel[[levelHere]] <- sort(unique(unlist(lapply(X = classesAndSubClasses, FUN = function(y){
			if(length(y) < levelHere) return(NULL)
			else return(y[[levelHere]])
		}))))
		labelsByLevel[[levelHere]] <- unlist(lapply(X = strsplit(x = classesByLevel[[levelHere]], split = "; "), FUN = tail, n=1))
	}
	
	## class counts
	countsByLevel <- list()
	for(levelHere in rev(seq_len(max(level)))){
		countsByLevel[[levelHere]] <- unlist(lapply(X = classesByLevel[[levelHere]], FUN = function(class){
			newSpectra <- ifelse(test = class %in% classifierClasses, yes = numberOfSpectra[[which(class == classifierClasses)]], no = 0)
			oldSpectra <- ifelse(test = levelHere < max(level), yes = sum(countsByLevel[[levelHere+1]][grepl(x = classesByLevel[[levelHere+1]], pattern = paste("^", class, sep = ""))]), no = 0)
			return(newSpectra + oldSpectra)
		}))
	}
	rootCount <- sum(countsByLevel[[1]])
	
	## coordinates
	colors <- rainbow(start = colorStart, alpha = colorAlpha, n = 1000)
	startDegreeByLevel <- list()
	spanDegreeByLevel <- list()
	colorByLevel <- list()
	
	for(levelHere in seq_len(max(level))){
		startDegreeByLevel[[levelHere]] <- list()
		spanDegreeByLevel[[levelHere]] <- list()
		colorByLevel[[levelHere]] <- list()
		
		classesToProcess <- classesByLevel[[levelHere]]
		precursorClasses <- NULL
		if(levelHere == 1)  precursorClasses <- ""
		else                precursorClasses <- classesByLevel[[levelHere-1]]
		
		for(precursorClassIdx in seq_along(precursorClasses)){
			precursorClass <- precursorClasses[[precursorClassIdx]]
			classesToProcessHereSelection <- grepl(x = classesToProcess, pattern = precursorClass)
			classesToProcessHere <- classesToProcess[classesToProcessHereSelection]
			startDegree <- ifelse(test = levelHere == 1, yes = 0, no = startDegreeByLevel[[levelHere-1]][[precursorClassIdx]])
			scalingFactor <- ifelse(test = levelHere == 1, yes = 1, no = countsByLevel[[levelHere-1]][[precursorClassIdx]] / sum(countsByLevel[[levelHere]][classesToProcessHereSelection])) ## ambiguous classes
			#startColor  <- ifelse(test = levelHere == 1, yes = 0, no = startDegreeByLevel[[levelHere-1]][[precursorClassIdx]])
			for(classToProcessHere in classesToProcessHere){
				classIdx <- which(classesByLevel[[levelHere]] == classToProcessHere)
				degreeSpan <- 360 * countsByLevel[[levelHere]][[classIdx]] / rootCount * scalingFactor
				startDegreeByLevel[[levelHere]][[classIdx]] <- startDegree
				spanDegreeByLevel [[levelHere]][[classIdx]] <- degreeSpan
				colorByLevel      [[levelHere]][[classIdx]] <- colors[[(floor(startDegree + degreeSpan / 2) / 360 * length(colors)) + ifelse(test = (floor(startDegree + degreeSpan / 2) / 360 * length(colors))==length(colors), yes = 0, no = 1) ]]
				startDegree <- startDegree + degreeSpan
			}
		}
	}
	
	thereIsNextLevelByLevel <- list()
	for(levelHere in seq_len(max(level))){
		thereIsNextLevelByLevel[[levelHere]] <- list()
		if(levelHere == max(level)){
			thereIsNextLevelByLevel[[levelHere]] <- rep(x = FALSE, times = length(classesByLevel[[levelHere]]))
		} else {
			for(classIdx in seq_along(classesByLevel[[levelHere]]))
				thereIsNextLevelByLevel[[levelHere]][[classIdx]] <- any(grepl(x = classesByLevel[[levelHere+1]], pattern = classesByLevel[[levelHere]][[classIdx]]))
		}
	}
	
	
	## If you use it in published research, please cite:
	##  Gu, Z. circlize implements and enhances circular visualization in R. Bioinformatics 2014.
	library("circlize")
	library("plotrix")
	
	degreeThresholdForDrawing <- 0.5
	#maxLevel <- max(which(unlist(lapply(X = spanDegreeByLevel, FUN = function(x){any(unlist(x) >= degreeThresholdForDrawing)}))))
	
	plotRadius <- 8
	plotCex1 <- 1
	plotCex2 <- 1
	
	plot(1, type="n", xlab="", ylab="", xlim=c(-plotRadius, plotRadius), ylim=c(-plotRadius, plotRadius), axes = FALSE)
	
	## circle segments
	tmp <- sapply(X = seq_along(classesByLevel), FUN = function(levelHere){
		sapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
			if(spanDegreeByLevel[[levelHere]][[classIdx]] < degreeThresholdForDrawing)  return()
			draw.sector(
				start.degree = startDegreeByLevel[[levelHere]][[classIdx]],
				end.degree = startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]],
				rou1 = levelHere - 1, rou2 = levelHere, center = c(0,0), clock.wise = FALSE, col = colorByLevel [[levelHere]][[classIdx]], border = "white"
			)
		})
	})
	## segment text
	minimumAngleToShowSegmentText <- 15
	tmp <- sapply(X = seq_along(classesByLevel), FUN = function(levelHere){
		sapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
			if(spanDegreeByLevel[[levelHere]][[classIdx]] < minimumAngleToShowSegmentText)  return()
			textTokens <- strwrap(x = labelsByLevel[[levelHere]][[classIdx]], width = max(nchar(strsplit(x = "Some text", split = " ")[[1]])))
			#firstOffset <- 1 / (length(textTokens) * 2)
			
			for(idx in seq_along(textTokens)){
				#offset <- firstOffset * (2 * idx - 1)
				middle <- (startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2)/360 * 2*pi
				isSwitched <- middle > pi/2 & middle < 3 * pi/2
				isSwitched <- middle > pi
				offset <-  ifelse(test = middle > pi, yes = (length(textTokens) - idx + 1) / (length(textTokens) + 1), no = idx / (length(textTokens) + 1))
				#if(isSwitched) textTokens[[idx]] <- rev(textTokens[[idx]])
				#label <- ifelse(test = isSwitched, yes = paste(strsplit(x = labelsByLevel[[levelHere]][[classIdx]], split = " ")[[1]], collapse = " "), no = labelsByLevel[[levelHere]][[classIdx]])
				#middle <- ifelse(test = middle < pi, yes = middle + pi, no = middle)
				arctext(
					x = textTokens[[idx]],
					center = c(0, 0), radius = levelHere - offset - 0.04,
					middle = middle,
					cex = plotCex1, stretch = 1,
					clockwise = !isSwitched
				)
			}
		})
	})
	
	## outer text
	levelMaxHere <- max(level)
	tmp <- sapply(X = seq_along(classesByLevel), FUN = function(levelHere){
		if(levelHere > levelMaxHere)  return()
		sapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
			if(thereIsNextLevelByLevel[[levelHere]][[classIdx]] | spanDegreeByLevel[[levelHere]][[classIdx]] >= minimumAngleToShowSegmentText)  return()
			#radius <- maxLevel + 0.2
			radius <- levelHere + 0.2
			angle <- (startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2)/360 * 2*pi
			x <- radius * sin(angle+pi/2)
			y <- radius * cos(angle+pi/2)
			srt <- startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2
			isSwitched <- srt > 90 & srt < 270
			if(isSwitched) adj <- c(1,0) else adj <- c(0,0.5)
			srt <- ifelse(test = isSwitched, yes = srt + 180, no = srt)
			text(
				x = x, y = -y, labels = labelsByLevel[[levelHere]][[classIdx]], adj = adj,
				srt=srt,
				cex = plotCex2
			)
		})
	})
}



# ---------- Export peak list as MAF ----------
f.export_maf <- function(peak_list, maf_filename) {
	# Preparations
	l <- nrow(peak_list)
	
	# These columns are defined by MetaboLights mzTab
	maf <- apply(X=data.frame(database_identifier=character(l),
							  chemical_formula=character(l),
							  smiles=character(l),
							  inchi=character(l),
							  metabolite_identification=character(l),
							  xcms_identifier=rownames(peak_list),
							  mass_to_charge=peak_list$mzmed,
							  fragmentation=character(l),
							  modifications=character(l),
							  charge=character(l),
							  retention_time=peak_list$rtmed,
							  taxid=character(l),
							  species=character(l),
							  database=character(l),
							  database_version=character(l),
							  reliability=character(l),
							  uri=character(l),
							  search_engine=character(l),
							  search_engine_score=character(l),
							  smallmolecule_abundance_sub=character(l),
							  smallmolecule_abundance_stdev_sub=character(l),
							  smallmolecule_abundance_std_error_sub=character(l),
							  peak_list,
							  stringsAsFactors=FALSE),
				 MARGIN=2, FUN=as.character)
	
	# Export MAF
	write.table(maf, file=maf_filename, row.names=FALSE, col.names=colnames(maf), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Update existing MAF with annotated compounds ----------
f.annotate_maf_compounds <- function(maf_input, maf_output, polarity, xcms_id, pol_mode, smiles, names) {
	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")
	
	maf_out$database_identifier <- as.character(maf_out$database_identifier)
	maf_out$metabolite_identification <- as.character(maf_out$metabolite_identification)
	
	# Annotate compounds
	for (i in c(1:length(xcms_id))) {
		if (pol_mode[i] == polarity) {
			j <- which(maf_out$xcms_identifier==xcms_id[i])
			if (length(j) > 0) {
				if (nchar(as.character(maf_out$database_identifier[j])) > 0) {
					maf_out$database_identifier[j] <- paste(c(smiles[i], maf_out$database_identifier[j]), collapse='|')
					maf_out$metabolite_identification[j] <- paste(c(names[i], maf_out$metabolite_identification[j]), collapse='|')
				} else {
					maf_out$database_identifier[j] <- paste(smiles[i], collapse='|')
					maf_out$metabolite_identification[j] <- paste(names[i], collapse='|')
				}
			}
		}
	}
	
	# Remove non-standard columns "ms_level", "primary_class"
	maf_out$ms_level <- NULL
	maf_out$primary_class <- NULL
	
	# Export MAF
	write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Update existing MAF with annotated compounds ----------
f.annotate_maf_classes <- function(maf_input, maf_output) {
	# Read CHEMONT ontology
	obo <- ontologyIndex::get_ontology(file=paste0("~/Desktop/Projekte/Habilitation/Mosses/ms-swath/data/ChemOnt_2_1.obo"), extract_tags="minimal")
	
	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")
	
	# Annotate classes
	for (i in c(1:length(maf_out$xcms_identifier))) {
		cl <- gsub(x=maf_out$primary_class[i], pattern='.*; ', replacement='')
		id <- as.character(obo$id[which(as.character(obo$name) %in% as.character(cl))])
		
		if (nchar(cl) > 1) {
			maf_out$database_identifier[i] <- paste(id, collapse='|')
			maf_out$metabolite_identification[i] <- paste(cl, collapse='|')
		}
	}
	
	# Export MAF
	write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Update existing MAF with annotated compounds ----------
f.annotate_maf_cached <- function(maf_input, maf_output, mzml_names, classes, classifiers) {
	# Annotate identified compounds
	max.mz.range <- 0.25
	max.rt.range <- 20
	
	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")
	
	# Read CHEMONT ontology
	obo <- ontologyIndex::get_ontology(file=paste0(mzml_dir, "../data/ChemOnt_2_1.obo"), extract_tags="minimal")
	
	# Get CHEMONTIDs for determined classes
	classes_id <- NULL
	for (i in 1:length(classes)) {
		classes_id <- c(classes_id, as.character(obo$id[which(as.character(obo$name) == gsub(x=classes[i], pattern=".*; ", replacement=""))]))
	}
	
	# Process annotated samples
	for (i in 1:length(mzml_names)) {
		for (j in 1:nrow(classifiers[[i]])) {
			# Compound info
			mz <- classifiers[[i]][j,"m.z"]
			rt <- classifiers[[i]][j,"RT"]
			cl <- gsub(x=classifiers[[i]][j,"Annotation..putative."], pattern=".*; ", replacement="")
			id <- as.character(obo$id[which(as.character(obo$name) == as.character(cl))])
			
			# Add info in MAF
			li <- which( (maf_out$mass_to_charge > (mz - max.mz.range)) & (maf_out$mass_to_charge < (mz + max.mz.range)) &
						 	(maf_out$retention_time > (rt - max.rt.range)) & (maf_out$retention_time < (rt + max.rt.range)) )
			if (length(li) == 0) {
				print(paste("No candidates in peak list found for",classifiers[[i]][j,"Label"],"in",classifiers[[i]][j,"Metabolite.name"]))
			} else {
				# Take the closest rt
				if (length(li) > 1) {
					li <- li[which.min(abs(maf_out[li,"retention_time"]-rt))]
				}
				# Write database_identifier
				#if (nchar(maf_out[li, "database_identifier"]) > 0)
				#	maf_out[li, "database_identifier"] <- paste(id, maf_out[li, "database_identifier"], sep=" | ")
				#else
					maf_out[li, "database_identifier"] <- id
				# Write metabolite_identification
				#if (nchar(maf_out[li, "metabolite_identification"]) > 0)
				#	maf_out[li, "metabolite_identification"] <- paste(cl, maf_out[li, "metabolite_identification"], sep=" | ")
				#else
					maf_out[li, "metabolite_identification"] <- cl
			}
		}
	}
	
	# Export MAF
	write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Update existing MAF with annotated compounds ----------
f.annotate_maf_each_pos <- function(maf_input, maf_output) {
	# Read CHEMONT ontology
	obo <- ontologyIndex::get_ontology(file=paste0(mzml_dir, "../data/ChemOnt_2_1.obo"), extract_tags="minimal")
	
	# Get CHEMONTIDs for determined classes
	classes_id <- NULL
	for (i in 1:length(classes_pos)) {
		classes_id <- c(classes_id, as.character(obo$id[which(as.character(obo$name) == gsub(x=classes_pos[i], pattern=".*; ", replacement=""))]))
	}
	
	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")
	
	# Fetch class for each spectrum
	for (i in maf_out$xcms_identifier) {
		#i<-"FT0383"
		print(i)
		# Create MSP entry
		RETENTIONTIME <- maf_out[which(maf_out$xcms_identifier==i), "rtmed"]
		msp_name <- paste0(i,".msp")
		msp_text <- NULL
		for (j in c(1,2)) {
			msp_text <- c(msp_text, paste("NAME:", i))
			msp_text <- c(msp_text, paste("AlignmentID:", j))
			msp_text <- c(msp_text, paste("RETENTIONTIME:", (RETENTIONTIME+1-j)))
			msp_text <- c(msp_text, paste("PRECURSORMZ:", maf_out[which(maf_out$xcms_identifier==i), "mzmed"]))
			msp_text <- c(msp_text, paste("METABOLITENAME:", i))
			if (feat_spectra_pos[[i]]@listData[[1]]@polarity == 1) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
			} else if (feat_spectra_pos[[i]]@listData[[1]]@polarity == 0) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
			} else {
				print("Warning: No adduction found.")
			}
			msp_text <- c(msp_text, paste("NumPeaks:", feat_spectra_pos[[i]]@listData[[1]]@peaksCount))
			PEAKS <- NULL
			for (k in 1:length(feat_spectra_pos[[i]]@listData[[1]]@mz)) {
				PEAKS <- c(PEAKS, paste(feat_spectra_pos[[i]]@listData[[1]]@mz[k], feat_spectra_pos[[i]]@listData[[1]]@intensity[k], sep="\t") )
			}
			msp_text <- c(msp_text, PEAKS)
			msp_text <- c(msp_text, "")
		}
		
		# Write MSP entry
		cat(msp_text, file=msp_name, sep="\n")
		
		# Get class for MSP entry
		classifier_object <- tryCatch(
			obj <- applyClassifierMs2(classifierFile = paste0("../data/", classifier_name_pos, ".RData"),
									  propertiesFile = paste0("../data/", classifier_name_pos, ".txt"),
									  fileMs1Path = NULL,
									  fileMs2Path = msp_name,
									  fileClasses = "../data/massbank_classes.txt",
									  minimumIntensityOfMaximalMS2peak = msms.intensity.threshold,
									  minimumProportionOfMS2peaks = 0.05,
									  mzDeviationAbsolute_grouping = mzabs,
									  mzDeviationInPPM_grouping = mzppm),
			error = function(c) NULL,
			warning = function(c) obj,
			message = function(c) obj)
		
		# Write class and id to maf object
		if (! is.null(classifier_object)) {
			cl <- gsub(x=unique(classifier_object$`Annotation (putative)`), pattern='.*; ', replacement='')
			id <- as.character(obo$id[which(as.character(obo$name) %in% as.character(cl))])
			
			maf_out$database_identifier[which(maf_out$xcms_identifier==i)] <- paste(id, collapse='|')
			maf_out$metabolite_identification[which(maf_out$xcms_identifier==i)] <- paste(cl, collapse='|')
		}
		
		# Remove MSP file
		file.remove(msp_name)
	}
	
	# Export MAF
	write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Update existing MAF with annotated compounds ----------
f.annotate_maf_each_parallel_pos <- function(maf_input, maf_output) {
	# Read CHEMONT ontology
	obo <- ontologyIndex::get_ontology(file=paste0(mzml_dir, "../data/ChemOnt_2_1.obo"), extract_tags="minimal")
	
	# Get CHEMONTIDs for determined classes
	classes_id <- NULL
	for (i in 1:length(classes_pos)) {
		classes_id <- c(classes_id, as.character(obo$id[which(as.character(obo$name) == gsub(x=classes_pos[i], pattern=".*; ", replacement=""))]))
	}
	
	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")
	
	# Fetch class for each spectrum
	maf_spectra <- list()
	maf_spectra <- foreach(i=c(as.character(maf_out$xcms_identifier))) %dopar% {
		# Create MSP entry
		print(i)
		RETENTIONTIME <- maf_out[which(maf_out$xcms_identifier==i), "rtmed"]
		msp_name <- paste0(i,".msp")
		msp_text <- NULL
		for (j in c(1,2)) {
			msp_text <- c(msp_text, paste("NAME:", i))
			msp_text <- c(msp_text, paste("AlignmentID:", j))
			msp_text <- c(msp_text, paste("RETENTIONTIME:", (RETENTIONTIME+1-j)))
			msp_text <- c(msp_text, paste("PRECURSORMZ:", maf_out[which(maf_out$xcms_identifier==i), "mzmed"]))
			msp_text <- c(msp_text, paste("METABOLITENAME:", i))
			if (feat_spectra_pos[[i]]@listData[[1]]@polarity == 1) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
			} else if (feat_spectra_pos[[i]]@listData[[1]]@polarity == 0) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
			} else {
				print("Warning: No adduction found.")
			}
			msp_text <- c(msp_text, paste("NumPeaks:", feat_spectra_pos[[i]]@listData[[1]]@peaksCount))
			PEAKS <- NULL
			for (k in 1:length(feat_spectra_pos[[i]]@listData[[1]]@mz)) {
				PEAKS <- c(PEAKS, paste(feat_spectra_pos[[i]]@listData[[1]]@mz[k], feat_spectra_pos[[i]]@listData[[1]]@intensity[k], sep="\t") )
			}
			msp_text <- c(msp_text, PEAKS)
			msp_text <- c(msp_text, "")
		}
		
		# Write MSP entry
		cat(msp_text, file=msp_name, sep="\n")
		
		# Get class for MSP entry
		classifier_object <- tryCatch(
			obj <- applyClassifierMs2(classifierFile = paste0("../data/", classifier_name_pos, ".RData"),
									  propertiesFile = paste0("../data/", classifier_name_pos, ".txt"),
									  fileMs1Path = NULL,
									  fileMs2Path = msp_name,
									  fileClasses = "../data/massbank_classes.txt",
									  minimumIntensityOfMaximalMS2peak = msms.intensity.threshold,
									  minimumProportionOfMS2peaks = 0.05,
									  mzDeviationAbsolute_grouping = mzabs,
									  mzDeviationInPPM_grouping = mzppm),
			error = function(c) NULL,
			warning = function(c) obj,
			message = function(c) obj)
		
		# Remove MSP file
		file.remove(msp_name)
		
		return(classifier_object)
	}
	
	for (i in c(1:length(maf_out$xcms_identifier))) {
		classifier_object <- maf_spectra[[i]]
		if (! is.null(classifier_object)) {
			cl <- gsub(x=unique(classifier_object$`Annotation (putative)`), pattern='.*; ', replacement='')
			id <- as.character(obo$id[which(as.character(obo$name) %in% as.character(cl))])
			
			maf_out$database_identifier[i] <- paste(id, collapse='|')
			maf_out$metabolite_identification[i] <- paste(cl, collapse='|')
		}
	}
	
	# Export MAF
	write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Update existing MAF with annotated compounds ----------
f.annotate_maf_each_neg <- function(maf_input, maf_output) {
	# Read CHEMONT ontology
	obo <- ontologyIndex::get_ontology(file=paste0(mzml_dir, "../data/ChemOnt_2_1.obo"), extract_tags="minimal")
	
	# Get CHEMONTIDs for determined classes
	classes_id <- NULL
	for (i in 1:length(classes_neg)) {
		classes_id <- c(classes_id, as.character(obo$id[which(as.character(obo$name) == gsub(x=classes_neg[i], pattern=".*; ", replacement=""))]))
	}
	
	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")
	
	# Fetch class for each spectrum
	for (i in maf_out$xcms_identifier) {
		#i<-"FT0383"
		print(i)
		# Create MSP entry
		RETENTIONTIME <- maf_out[which(maf_out$xcms_identifier==i), "rtmed"]
		msp_name <- paste0(i,".msp")
		msp_text <- NULL
		for (j in c(1,2)) {
			msp_text <- c(msp_text, paste("NAME:", i))
			msp_text <- c(msp_text, paste("AlignmentID:", j))
			msp_text <- c(msp_text, paste("RETENTIONTIME:", (RETENTIONTIME+1-j)))
			msp_text <- c(msp_text, paste("PRECURSORMZ:", maf_out[which(maf_out$xcms_identifier==i), "mzmed"]))
			msp_text <- c(msp_text, paste("METABOLITENAME:", i))
			if (feat_spectra_neg[[i]]@listData[[1]]@polarity == 1) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
			} else if (feat_spectra_neg[[i]]@listData[[1]]@polarity == 0) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
			} else {
				print("Warning: No adduction found.")
			}
			msp_text <- c(msp_text, paste("NumPeaks:", feat_spectra_neg[[i]]@listData[[1]]@peaksCount))
			PEAKS <- NULL
			for (k in 1:length(feat_spectra_neg[[i]]@listData[[1]]@mz)) {
				PEAKS <- c(PEAKS, paste(feat_spectra_neg[[i]]@listData[[1]]@mz[k], feat_spectra_neg[[i]]@listData[[1]]@intensity[k], sep="\t") )
			}
			msp_text <- c(msp_text, PEAKS)
			msp_text <- c(msp_text, "")
		}
		
		# Write MSP entry
		cat(msp_text, file=msp_name, sep="\n")
		
		# Get class for MSP entry
		classifier_object <- tryCatch(
			obj <- applyClassifierMs2(classifierFile = paste0("../data/", classifier_name_neg, ".RData"),
									  propertiesFile = paste0("../data/", classifier_name_neg, ".txt"),
									  fileMs1Path = NULL,
									  fileMs2Path = msp_name,
									  fileClasses = "../data/massbank_classes.txt",
									  minimumIntensityOfMaximalMS2peak = msms.intensity.threshold,
									  minimumProportionOfMS2peaks = 0.05,
									  mzDeviationAbsolute_grouping = mzabs,
									  mzDeviationInPPM_grouping = mzppm),
			error = function(c) NULL,
			warning = function(c) obj,
			message = function(c) obj)
		
		# Write class and id to maf object
		if (! is.null(classifier_object)) {
			cl <- gsub(x=unique(classifier_object$`Annotation (putative)`), pattern='.*; ', replacement='')
			id <- as.character(obo$id[which(as.character(obo$name) %in% as.character(cl))])
			
			maf_out$database_identifier[which(maf_out$xcms_identifier==i)] <- paste(id, collapse='|')
			maf_out$metabolite_identification[which(maf_out$xcms_identifier==i)] <- paste(cl, collapse='|')
		}
		
		# Remove MSP file
		file.remove(msp_name)
	}
	
	# Export MAF
	write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Update existing MAF with annotated compounds ----------
f.annotate_maf_each_parallel_neg <- function(maf_input, maf_output) {
	# Read CHEMONT ontology
	obo <- ontologyIndex::get_ontology(file=paste0(mzml_dir, "../data/ChemOnt_2_1.obo"), extract_tags="minimal")
	
	# Get CHEMONTIDs for determined classes
	classes_id <- NULL
	for (i in 1:length(classes_neg)) {
		classes_id <- c(classes_id, as.character(obo$id[which(as.character(obo$name) == gsub(x=classes_neg[i], pattern=".*; ", replacement=""))]))
	}
	
	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")
	
	# Fetch class for each spectrum
	maf_spectra <- list()
	maf_spectra <- foreach(i=c(as.character(maf_out$xcms_identifier))) %dopar% {
		# Create MSP entry
		print(i)
		RETENTIONTIME <- maf_out[which(maf_out$xcms_identifier==i), "rtmed"]
		msp_name <- paste0(i,".msp")
		msp_text <- NULL
		for (j in c(1,2)) {
			msp_text <- c(msp_text, paste("NAME:", i))
			msp_text <- c(msp_text, paste("AlignmentID:", j))
			msp_text <- c(msp_text, paste("RETENTIONTIME:", (RETENTIONTIME+1-j)))
			msp_text <- c(msp_text, paste("PRECURSORMZ:", maf_out[which(maf_out$xcms_identifier==i), "mzmed"]))
			msp_text <- c(msp_text, paste("METABOLITENAME:", i))
			if (feat_spectra_neg[[i]]@listData[[1]]@polarity == 1) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
			} else if (feat_spectra_neg[[i]]@listData[[1]]@polarity == 0) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
			} else {
				print("Warning: No adduction found.")
			}
			msp_text <- c(msp_text, paste("NumPeaks:", feat_spectra_neg[[i]]@listData[[1]]@peaksCount))
			PEAKS <- NULL
			for (k in 1:length(feat_spectra_neg[[i]]@listData[[1]]@mz)) {
				PEAKS <- c(PEAKS, paste(feat_spectra_neg[[i]]@listData[[1]]@mz[k], feat_spectra_neg[[i]]@listData[[1]]@intensity[k], sep="\t") )
			}
			msp_text <- c(msp_text, PEAKS)
			msp_text <- c(msp_text, "")
		}
		
		# Write MSP entry
		cat(msp_text, file=msp_name, sep="\n")
		
		# Get class for MSP entry
		classifier_object <- tryCatch(
			obj <- applyClassifierMs2(classifierFile = paste0("../data/", classifier_name_neg, ".RData"),
									  propertiesFile = paste0("../data/", classifier_name_neg, ".txt"),
									  fileMs1Path = NULL,
									  fileMs2Path = msp_name,
									  fileClasses = "../data/massbank_classes.txt",
									  minimumIntensityOfMaximalMS2peak = msms.intensity.threshold,
									  minimumProportionOfMS2peaks = 0.05,
									  mzDeviationAbsolute_grouping = mzabs,
									  mzDeviationInPPM_grouping = mzppm),
			error = function(c) NULL,
			warning = function(c) obj,
			message = function(c) obj)
		
		# Remove MSP file
		file.remove(msp_name)
		
		return(classifier_object)
	}
	
	for (i in c(1:length(maf_out$xcms_identifier))) {
		classifier_object <- maf_spectra[[i]]
		if (! is.null(classifier_object)) {
			cl <- gsub(x=unique(classifier_object$`Annotation (putative)`), pattern='.*; ', replacement='')
			id <- as.character(obo$id[which(as.character(obo$name) %in% as.character(cl))])
			
			maf_out$database_identifier[i] <- paste(id, collapse='|')
			maf_out$metabolite_identification[i] <- paste(cl, collapse='|')
		}
	}
	
	# Export MAF
	write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Classify only selected features in pos mode ----------
f.classify_selected_pos <- function(selected_features, maf_input) {
	# Read CHEMONT ontology
	obo <- ontologyIndex::get_ontology(file=paste0(mzml_dir, "../data/ChemOnt_2_1.obo"), extract_tags="minimal")
	
	# Get CHEMONTIDs for determined classes
	classes_id <- NULL
	for (i in 1:length(classes_pos)) {
		classes_id <- c(classes_id, as.character(obo$id[which(as.character(obo$name) == gsub(x=classes_pos[i], pattern=".*; ", replacement=""))]))
	}
	
	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")
	
	# Fetch class for each spectrum
	maf_spectra <- list()
	maf_spectra <- foreach(i=c(as.character(selected_features))) %dopar% {
		# Create MSP entry
		print(i)
		RETENTIONTIME <- maf_out[which(maf_out$xcms_identifier==i), "rtmed"]
		msp_name <- paste0(i,".msp")
		msp_text <- NULL
		for (j in c(1,2)) {
			msp_text <- c(msp_text, paste("NAME:", i))
			msp_text <- c(msp_text, paste("AlignmentID:", j))
			msp_text <- c(msp_text, paste("RETENTIONTIME:", (RETENTIONTIME+1-j)))
			msp_text <- c(msp_text, paste("PRECURSORMZ:", maf_out[which(maf_out$xcms_identifier==i), "mzmed"]))
			msp_text <- c(msp_text, paste("METABOLITENAME:", i))
			if (feat_spectra_pos[[i]]@listData[[1]]@polarity == 1) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
			} else if (feat_spectra_pos[[i]]@listData[[1]]@polarity == 0) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
			} else {
				print("Warning: No adduction found.")
			}
			msp_text <- c(msp_text, paste("NumPeaks:", feat_spectra_pos[[i]]@listData[[1]]@peaksCount))
			PEAKS <- NULL
			for (k in 1:length(feat_spectra_pos[[i]]@listData[[1]]@mz)) {
				PEAKS <- c(PEAKS, paste(feat_spectra_pos[[i]]@listData[[1]]@mz[k], feat_spectra_pos[[i]]@listData[[1]]@intensity[k], sep="\t") )
			}
			msp_text <- c(msp_text, PEAKS)
			msp_text <- c(msp_text, "")
		}
		
		# Write MSP entry
		cat(msp_text, file=msp_name, sep="\n")
		
		# Get class for MSP entry
		classifier_object <- tryCatch(
			obj <- applyClassifierMs2(classifierFile = paste0("../data/", classifier_name_pos, ".RData"),
									  propertiesFile = paste0("../data/", classifier_name_pos, ".txt"),
									  fileMs1Path = NULL,
									  fileMs2Path = msp_name,
									  fileClasses = "../data/massbank_classes.txt",
									  minimumIntensityOfMaximalMS2peak = msms.intensity.threshold,
									  minimumProportionOfMS2peaks = 0.05,
									  mzDeviationAbsolute_grouping = mzabs,
									  mzDeviationInPPM_grouping = mzppm),
			error = function(c) NULL,
			warning = function(c) obj,
			message = function(c) obj)
		
		# Remove MSP file
		file.remove(msp_name)
		
		return(classifier_object)
	}
	
	# Diversity of selected classes
	selected_classes_pos <- NULL
	for (i in c(1:length(maf_spectra))) {
		classifier_object <- maf_spectra[[i]]
		if (! is.null(classifier_object)) {
			for (j in c(1:length(unique(classifier_object$`Annotation (putative)`)))) {
				if (gsub(x=unique(classifier_object$`Annotation (putative)`)[j], pattern='.*; ', replacement='') != "") {
					selected_classes_pos <- c(selected_classes_pos, unique(classifier_object$`Annotation (putative)`)[j])
				}
			}
		}
	}
	div_selected_classes_pos <- as.data.frame(table(selected_classes_pos))
	rownames(div_selected_classes_pos) <- div_selected_classes_pos$selected_classes_pos
	div_selected_classes_pos$selected_classes_pos <- NULL
	
	# div_selected_classes_samples
	div_selected_classes_samples_pos <- data.frame(selected_classes_pos=div_selected_classes_pos$selected_classes_pos)
	for (i in mzml_names_pos) div_selected_classes_samples_pos <- cbind(div_selected_classes_samples_pos, data.frame(rep(0, nrow(div_selected_classes_samples_pos))))
	rownames(div_selected_classes_samples_pos) <- div_selected_classes_samples_pos$selected_classes_pos
	div_selected_classes_samples_pos$selected_classes_pos <- NULL
	colnames(div_selected_classes_samples_pos) <- mzml_names_pos
	
	for (i in mzml_names_pos) {
		for (j in names(which(bina_list_pos[rownames(bina_list_pos)==i, colnames(bina_list_pos) %in% as.character(selected_features)] != 0))) {
			classifier_object <- maf_spectra[[which(selected_features==j)]]
			for (k in c(1:length(classifier_object$`Annotation (putative)`))) {
				if (! is.null(classifier_object$`Annotation (putative)`[k])) {
					if (gsub(x=classifier_object$`Annotation (putative)`[k], pattern='.*; ', replacement='') != "") {
						div_selected_classes_samples_pos[rownames(div_selected_classes_samples_pos)==classifier_object$`Annotation (putative)`[k], colnames(div_selected_classes_samples_pos)==i] <- div_selected_classes_samples_pos[rownames(div_selected_classes_samples_pos)==classifier_object$`Annotation (putative)`[k], colnames(div_selected_classes_samples_pos)==i] + 1
					}
				}
			}
		}
	}
	
	# Diversity of classes per species
	div_selected_classes_species_pos <- NULL
	for (i in species_names) {
		obj <- div_selected_classes_samples_pos[, which(mzml_names_pos %in% mzml_names_pos[which(species==i)]) ]
		obj[is.na(obj)] <- 0
		obj <- apply(X=obj, MARGIN=1, FUN=function(x) { sum(x) })
		obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
		obj[which(obj$frequency==0), "frequency"] <- 0
		if (is.null(div_selected_classes_species_pos)) {
			div_selected_classes_species_pos <- obj
		} else {
			div_selected_classes_species_pos <- merge(div_selected_classes_species_pos, obj, by="classes", all.x=TRUE, all.y=TRUE)
		}
	}
	rownames(div_selected_classes_species_pos) <- div_selected_classes_species_pos$classes
	div_selected_classes_species_pos <- div_selected_classes_species_pos[, -which(colnames(div_selected_classes_species_pos)=="classes")]
	colnames(div_selected_classes_species_pos) <- species_names
	
	# Save global variables
	selected_classes_pos <<- selected_classes_pos
	div_selected_classes_pos <<- div_selected_classes_pos
	div_selected_classes_samples_pos <<- div_selected_classes_samples_pos
	div_selected_classes_species_pos <<- div_selected_classes_species_pos
	
	# Return classification result
	return(maf_spectra)
}



# ---------- Classify only selected features in neg mode ----------
f.classify_selected_neg <- function(selected_features, maf_input) {
	# Read CHEMONT ontology
	obo <- ontologyIndex::get_ontology(file=paste0(mzml_dir, "../data/ChemOnt_2_1.obo"), extract_tags="minimal")
	
	# Get CHEMONTIDs for determined classes
	classes_id <- NULL
	for (i in 1:length(classes_neg)) {
		classes_id <- c(classes_id, as.character(obo$id[which(as.character(obo$name) == gsub(x=classes_neg[i], pattern=".*; ", replacement=""))]))
	}
	
	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")
	
	# Fetch class for each spectrum
	maf_spectra <- list()
	maf_spectra <- foreach(i=c(as.character(selected_features))) %dopar% {
		# Create MSP entry
		print(i)
		RETENTIONTIME <- maf_out[which(maf_out$xcms_identifier==i), "rtmed"]
		msp_name <- paste0(i,".msp")
		msp_text <- NULL
		for (j in c(1,2)) {
			msp_text <- c(msp_text, paste("NAME:", i))
			msp_text <- c(msp_text, paste("AlignmentID:", j))
			msp_text <- c(msp_text, paste("RETENTIONTIME:", (RETENTIONTIME+1-j)))
			msp_text <- c(msp_text, paste("PRECURSORMZ:", maf_out[which(maf_out$xcms_identifier==i), "mzmed"]))
			msp_text <- c(msp_text, paste("METABOLITENAME:", i))
			if (feat_spectra_neg[[i]]@listData[[1]]@polarity == 1) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
			} else if (feat_spectra_neg[[i]]@listData[[1]]@polarity == 0) {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
			} else {
				print("Warning: No adduction found.")
			}
			msp_text <- c(msp_text, paste("NumPeaks:", feat_spectra_neg[[i]]@listData[[1]]@peaksCount))
			PEAKS <- NULL
			for (k in 1:length(feat_spectra_neg[[i]]@listData[[1]]@mz)) {
				PEAKS <- c(PEAKS, paste(feat_spectra_neg[[i]]@listData[[1]]@mz[k], feat_spectra_neg[[i]]@listData[[1]]@intensity[k], sep="\t") )
			}
			msp_text <- c(msp_text, PEAKS)
			msp_text <- c(msp_text, "")
		}
		
		# Write MSP entry
		cat(msp_text, file=msp_name, sep="\n")
		
		# Get class for MSP entry
		classifier_object <- tryCatch(
			obj <- applyClassifierMs2(classifierFile = paste0("../data/", classifier_name_neg, ".RData"),
									  propertiesFile = paste0("../data/", classifier_name_neg, ".txt"),
									  fileMs1Path = NULL,
									  fileMs2Path = msp_name,
									  fileClasses = "../data/massbank_classes.txt",
									  minimumIntensityOfMaximalMS2peak = msms.intensity.threshold,
									  minimumProportionOfMS2peaks = 0.05,
									  mzDeviationAbsolute_grouping = mzabs,
									  mzDeviationInPPM_grouping = mzppm),
			error = function(c) NULL,
			warning = function(c) obj,
			message = function(c) obj)
		
		# Remove MSP file
		file.remove(msp_name)
		
		return(classifier_object)
	}
	
	# Diversity of selected classes
	selected_classes_neg <- NULL
	for (i in c(1:length(maf_spectra))) {
		classifier_object <- maf_spectra[[i]]
		if (! is.null(classifier_object)) {
			for (j in c(1:length(unique(classifier_object$`Annotation (putative)`)))) {
				if (gsub(x=unique(classifier_object$`Annotation (putative)`)[j], pattern='.*; ', replacement='') != "") {
					selected_classes_neg <- c(selected_classes_neg, unique(classifier_object$`Annotation (putative)`)[j])
				}
			}
		}
	}
	div_selected_classes_neg <- as.data.frame(table(selected_classes_neg))
	rownames(div_selected_classes_neg) <- div_selected_classes_neg$selected_classes_neg
	div_selected_classes_neg$selected_classes_neg <- NULL
	
	# div_selected_classes_samples
	div_selected_classes_samples_neg <- data.frame(selected_classes_neg=div_selected_classes_neg$selected_classes_neg)
	for (i in mzml_names_neg) div_selected_classes_samples_neg <- cbind(div_selected_classes_samples_neg, data.frame(rep(0, nrow(div_selected_classes_samples_neg))))
	rownames(div_selected_classes_samples_neg) <- div_selected_classes_samples_neg$selected_classes_neg
	div_selected_classes_samples_neg$selected_classes_neg <- NULL
	colnames(div_selected_classes_samples_neg) <- mzml_names_neg
	
	for (i in mzml_names_neg) {
		for (j in names(which(bina_list_neg[rownames(bina_list_neg)==i, colnames(bina_list_neg) %in% as.character(selected_features)] != 0))) {
			classifier_object <- maf_spectra[[which(selected_features==j)]]
			for (k in c(1:length(classifier_object$`Annotation (putative)`))) {
				if (! is.null(classifier_object$`Annotation (putative)`[k])) {
					if (gsub(x=classifier_object$`Annotation (putative)`[k], pattern='.*; ', replacement='') != "") {
						div_selected_classes_samples_neg[rownames(div_selected_classes_samples_neg)==classifier_object$`Annotation (putative)`[k], colnames(div_selected_classes_samples_neg)==i] <- div_selected_classes_samples_neg[rownames(div_selected_classes_samples_neg)==classifier_object$`Annotation (putative)`[k], colnames(div_selected_classes_samples_neg)==i] + 1
					}
				}
			}
		}
	}
	
	# Diversity of classes per species
	div_selected_classes_species_neg <- NULL
	for (i in species_names) {
		obj <- div_selected_classes_samples_neg[, which(mzml_names_neg %in% mzml_names_neg[which(species==i)]) ]
		obj[is.na(obj)] <- 0
		obj <- apply(X=obj, MARGIN=1, FUN=function(x) { sum(x) })
		obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
		obj[which(obj$frequency==0), "frequency"] <- 0
		if (is.null(div_selected_classes_species_neg)) {
			div_selected_classes_species_neg <- obj
		} else {
			div_selected_classes_species_neg <- merge(div_selected_classes_species_neg, obj, by="classes", all.x=TRUE, all.y=TRUE)
		}
	}
	rownames(div_selected_classes_species_neg) <- div_selected_classes_species_neg$classes
	div_selected_classes_species_neg <- div_selected_classes_species_neg[, -which(colnames(div_selected_classes_species_neg)=="classes")]
	colnames(div_selected_classes_species_neg) <- species_names
	
	# Save global variables
	selected_classes_neg <<- selected_classes_neg
	div_selected_classes_neg <<- div_selected_classes_neg
	div_selected_classes_samples_neg <<- div_selected_classes_samples_neg
	div_selected_classes_species_neg <<- div_selected_classes_species_neg
	
	# Return classification result
	return(maf_spectra)
}



# ---------- Annotate MetFamily Classification Performances in pos mode ----------
f.annotate_metfamily_classification_performances_pos <- function() {
	# Read CHEMONT ontology
	obo <- ontologyIndex::get_ontology(file=paste0(mzml_dir, "../data/ChemOnt_2_1.obo"), extract_tags="minimal")
	
	# Read MetFamily performances in pos mode
	metfamily_perf_pos <- read.table(file="data/ms2_pos_classified_performance.csv", quote="\"", sep=";", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	
	# Get CHEMONTIDs for determined classes
	classes_id <- NULL
	for (i in 1:length(metfamily_perf_pos$compound_class)) {
		classes_id <- c(classes_id, gsub(x=as.character(obo$id[which(as.character(obo$name) == gsub(x=metfamily_perf_pos$compound_class[i], pattern=".*; ", replacement=""))]),
										 pattern='CHEMONTID:', replacement=''))
	}
	
	metfamily_perf_pos$compound_class <- gsub(x=metfamily_perf_pos$compound_class, pattern=".*; ", replacement="")
	metfamily_perf_pos$chemont.id <- classes_id
	
	# Write table to vignettes
	write.table(metfamily_perf_pos, file="vignettes/ms2_pos_classified_performance.csv", row.names=FALSE, col.names=colnames(metfamily_perf_pos), quote=TRUE, sep=";", na="\"\"")
}



# ---------- Annotate MetFamily Classification Performances in neg mode ----------
f.annotate_metfamily_classification_performances_neg <- function() {
	# Read CHEMONT ontology
	obo <- ontologyIndex::get_ontology(file=paste0(mzml_dir, "../data/ChemOnt_2_1.obo"), extract_tags="minimal")
	
	# Read MetFamily performances in neg mode
	metfamily_perf_neg <- read.table(file="data/ms2_neg_classified_performance.csv", quote="\"", sep=";", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	
	# Get CHEMONTIDs for determined classes
	classes_id <- NULL
	for (i in 1:length(metfamily_perf_neg$compound_class)) {
		classes_id <- c(classes_id, gsub(x=as.character(obo$id[which(as.character(obo$name) == gsub(x=metfamily_perf_neg$compound_class[i], pattern=".*; ", replacement=""))]),
										 pattern='CHEMONTID:', replacement=''))
	}
	
	metfamily_perf_neg$compound_class <- gsub(x=metfamily_perf_neg$compound_class, pattern=".*; ", replacement="")
	metfamily_perf_neg$chemont.id <- classes_id
	
	# Write table to vignettes
	write.table(metfamily_perf_neg, file="vignettes/ms2_neg_classified_performance.csv", row.names=FALSE, col.names=colnames(metfamily_perf_neg), quote=TRUE, sep=";", na="\"\"")
}



# ---------- Query Online COCONUT ----------
f.query_coconut <- function(smiles) {
	#library(httr)
	#library(jsonlite)
	#smiles="O=C1OC(C(O)=C1O)CO"
	
	# Construct URL
	url <- paste0("https://coconut.naturalproducts.net/api/search/exact-structure?type=inchi&smiles=", smiles)
	
	# Make query to COCONUT
	query <- httr::GET(url=url)
	
	# Export data frame
	if (query$status_code >= 400) {
		httr::stop_for_status(x=query, task=paste0("Error! Query to COCONUT failed with error code #", query$status_code))
		return(data.frame())
	} else {
		result <- httr::content(x=query, as="text")
		object <- fromJSON(txt=result, flatten=TRUE)
		
		if (object$count > 0) {
			return(object$naturalProducts)
		} else {
			return(data.frame())
		}
	}
}



# ---------- Query Online LOTUS ----------
f.query_lotus <- function(smiles) {
	#library(httr)
	#library(jsonlite)
	#smiles="O=C1OC(C(O)=C1O)CO"
	
	# Construct URL
	url <- paste0("https://lotus.naturalproducts.net/api/search/exact-structure?type=inchi&smiles=", smiles)
	
	# Make query to LOTUS
	query <- httr::GET(url=url)
	
	# Export data frame
	if (query$status_code >= 400) {
		httr::stop_for_status(x=query, task=paste0("Error! Query to COCONUT failed with error code #", query$status_code))
		return(data.frame())
	} else {
		result <- httr::content(x=query, as="text")
		object <- fromJSON(txt=result, flatten=TRUE)
		
		if (object$count > 0) {
			return(object$naturalProducts)
		} else {
			return(data.frame())
		}
	}
}



# ---------- Query Online LOTUS ----------
f.make_classes_at_chemont_level <- function(try_superclass_level=3, div_classes_samples) {
	try_classes_samples_names <- NULL
	for (i in c(1:try_superclass_level)) {
		try_classes_samples_names <- c(try_classes_samples_names, lapply(X=strsplit(rownames(div_classes_samples), '; '), FUN=function(x) { gsub(x=paste(x[1:i],sep='',collapse='; '),pattern='; NA',replacement='') }))
	}
	try_classes_samples <- data.frame()
	for (i in c(1:ncol(div_classes_samples))) try_classes_samples <- rbind(try_classes_samples, rep(0, length(unique(try_classes_samples_names))))
	try_classes_samples <- t(try_classes_samples)
	colnames(try_classes_samples) <- colnames(div_classes_samples)
	rownames(try_classes_samples) <- unique(try_classes_samples_names)
	for (i in rownames(try_classes_samples)) {
		for (j in c(1:ncol(div_classes_samples))) {
			try_classes_samples[rownames(try_classes_samples)==i, j] <- sum(div_classes_samples[grep(x=rownames(div_classes_samples), pattern=i), j])
		}
	}
	
	# Diversity of superclasses
	try_classes <- try_classes_samples
	try_classes[is.na(try_classes)] <- 0
	try_classes <- apply(X=try_classes, MARGIN=1, FUN=function(x) { sum(x) })
	try_classes <- data.frame(row.names=names(try_classes), frequency=as.numeric(try_classes))
	
	# Imputation of NA with zeros
	try_classes[is.na(try_classes)] <- 0
	try_classes_samples[is.na(try_classes_samples)] <- 0
	
	# Classification list for statistics
	try_class_list <- as.data.frame(t(try_classes_samples))
	try_class_list[is.na(try_class_list)] <- 0
	
	# Return
	return(try_class_list)
}


