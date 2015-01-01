## script to load in latest EWAS results
## mortality, diabetes, preterm birth, and quantitative phenotypes
library('rmeta')
combine_estimates <- function(validatedData) {
	varnames <- unique(validatedData$varname)
	combinedData <- data.frame()
	for(ii in 1:length(varnames)) {
		dat <- subset(validatedData, varname == varnames[ii])
		summ <- dat$estimate[!is.na(dat$estimate)]
		ses <- dat$standard_error[!is.na(dat$estimate)]
		effect <- summ
		pval <- dat$pvalue;
		se <- ses
		if(length(summ) > 1) {
			m <- meta.summaries(summ, ses, method='random')
			effect <- m$summary
			se <- m$se.summary
			pval <- pnorm(abs(effect/se), lower.tail=F)
		} 
		combinedData <- rbind(combinedData, data.frame(varname=varnames[ii], effect_size=effect, pvalue_combined=pval))
	}
	return(combinedData)
}

# 1.) mortality
load('../../result/ewas_results/ewas_mortality.Rdata')
ewasData.mortality <- allData

# 2.) diabetes
load('../../result/ewas_results/ewas_diabetes.Rdata')
sigTab <- table(subset(allData, pvalue <= 0.05)$varname)
allData$validated <- 0 
allData[allData$varname %in% names(sigTab)[sigTab>=2], 'validated'] <- 1
combined.diabetes <- combine_estimates(allData)
ewasData.diabetes <- combined.diabetes
ewasData.diabetes$validated <- 0
ewasData.diabetes$validated[ewasData.diabetes$varname %in% allData$varname[allData$validated == 1]] <- 1
ewasData.diabetes$validated[ewasData.diabetes$varname == 'LBXHBS'] <- 0

# 3.) preterm birth
load('../../result/ewas_results/ewas_preterm.Rdata')
ewasData.preterm <- allData
ewasData.preterm$validated <- 0
ewasData.preterm$validated[ewasData.preterm$varname == 'URXBPH'] <- 1

# 4.) quantitative phenotypes
load('../../result/ewas_results/ewas_clinical_phenome.Rdata')
ewasData.quant <- metaEwas

