


# ---------- Shannon Diversity ----------
#' Function to calculate the Shannon diversity measure
#'
#' @param p vector with response variables of one sample
#' @export
#' @import vegan MASS
#' @examples
#' shannon.diversity(p=c(4,8))
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
#' Function to calculate the Menhinick diversity measure
#'
#' @param p vector with response variables of one sample
#' @export
#' @import vegan MASS
#' @examples
#' menhinick.diversity(p=c(4,8))
menhinick.diversity <- function(p) {
	# Based on: http://www.coastalwiki.org/wiki/Measurements_of_biodiversity#Species_richness_indices
	D_Mn <- length(p) / sqrt(vegan::specnumber(p))
}



# ---------- Tukey-Test ----------
#' Function to perform the Tukey post-hoc HSD test on response and model terms
#'
#' @param response vector with response variables
#' @param term vector with factorized terms
#' @export
#' @import stats multcomp
###_ @examples
###_ tukey.test(response=model_div$unique, term=as.factor(mzml_pheno_samples))
tukey.test <- function(response, term) {
	model_anova <- stats::aov(formula(response ~ term))
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



