## Chirag Patel
## 1/14/14
## correlation engine 
library(survey)
library(getopt)
options(survey.lonely.psu="adjust")

ADJUST_VARS <- c()
### main data file
#load('~/EWAS/LD/corrdata6/bigtable_correlation_nhanes.Rdata')
#load('~/projects/exposome_corrs/data/corrdata6/bigtable_correlation_nhanes.Rdata')
#load('/home/cjp/bigtable_correlation_nhanes.Rdata')

load('../data/corrdata6/bigtable_correlation_nhanes.Rdata')

spec = matrix(c('outpath', 'o', 1, 'character',
				'series', 's', 1, 'character',
				'variable', 'v', 1, 'character',
				'permute', 'p', 1, 'integer',
				'parametric', 'c', 1, 'integer',
				'strata', 'l', 1, 'character'
				),
				nrow=6, byrow=TRUE);
opt <- getopt(spec)
## parameters: variable, survey year, data source script, outfolder, sql for table, and permute

DATASCRIPT <- NULL
STRATA <- NULL
if(is.null(opt$series)) {SERIESNAME <- '2003-2004';} else {SERIESNAME <- opt$series;}
if(is.null(opt$variable)) {variable <- 'URXUCD';} else {variable <- opt$variable;}
if(is.null(opt$permute)) {PERMUTE <- 0;} else {PERMUTE <- opt$permute;}
if(is.null(opt$parametric)) {PARAMETRIC <- TRUE;} else {PARAMETRIC <- opt$parametric == 1;}
if(is.null(opt$outpath)) {OUTFOLDER <- '.';} else {OUTFOLDER <- opt$outpath;}
if(!is.null(opt$strata)) {STRATA <- opt$strata;} 
#if(!is.null(opt$script)) {DATASCRIPT <- opt$script;} 


outfilename <- sprintf('%s_%s.Rdata', variable, SERIESNAME)
if(PERMUTE > 0) {
	outfilename <- sprintf('%s_%s_permute.Rdata', variable, SERIESNAME)
}

filename <- file.path(OUTFOLDER, outfilename)
cat(sprintf('variable/series:%s/%s\n', variable, SERIESNAME))
cat(sprintf('permuted run:%i\n', PERMUTE))
cat(sprintf('output filename:%s\n', filename))
cat(sprintf('pearson(1)/spearman(0)):%i\n', PARAMETRIC))
#cat(sprintf('data script:%s\n', DATASCRIPT))
cat(sprintf('adjustment variables (partial correlation):%s\n', paste(ADJUST_VARS, collapse=",")))


svycorr <- function(var1,var2,dsn) {
	##weighted pearson corr
	form <- as.formula(paste("~", var1, "+", var2, sep=""))
	covr <- svyvar(form, dsn,na.rm=T)
	corrl <- NA
	model<- NA
	if(ncol(covr) == 2 & nrow(covr) == 2) {
		corrl <- covr[1,2] / (sqrt(covr[1,1]) * sqrt(covr[2,2]))
	}
	N <- sum(!is.na(dsn$variables[,var1]) & !is.na(dsn$variables[,var2]))
	return(list(cor=corrl, covar=covr, N=N))
}

pearsoncorr <- function(var1, var2, dat) {
	# unweighted pearson corr
	corrl <- cor( dat[, var1],  dat[, var2], use='pairwise.complete.obs', method = "pearson")
	N <- sum(!is.na(dat[,var1]) & !is.na(dat[,var2]))
	#corrl <- corr.test(dat[, var1], dat[,var2] , method='pearson', adjust='none')
	return(list(cor=corrl, covar=NA, N=N))
	#return(list(cor=corrl$r, covar=NA, N=corrl$n,p=corrl$p, t=corrl$t))
}

spearmancorr <- function(var1, var2, dat) {
	corrl <- cor( dat[, var1],  dat[, var2], use='pairwise.complete.obs', method = "spearman")
	N <- sum(!is.na(dat[,var1]) & !is.na(dat[,var2]))
	#corrl <- corr.test(dat[, var1], dat[,var2] , method='spearman', adjust='none')
	return(list(cor=corrl, covar=NA, N=N))
	#return(list(cor=corrl$r, covar=NA, N=corrl$n,p=corrl$p, t=corrl$t))
}


regressioncorr <- function(var1, var2, dsn) {
	frm <- as.formula(sprintf('%s ~ %s', var1, var2))
	mod <- svyglm(frm, dsn)
	summ <- coef(summary(mod))
	frmvar <- as.formula(sprintf('~ %s + %s', var1, var2))
	vv <- as.matrix(svyvar(frmvar, dsn, na.rm=T))
	sd1 <- sqrt(vv[1,1])
	sd2 <- sqrt(vv[2,2])	
	corrl <- (summ[2, 1]*sd2)/sd1
	pval <- summ[2, 4]
	return(list(cor=corrl, pvalue=pval))
}

surveyDesign <- function(dat) {
	dsn <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weight=~WTMEC2YR, data=subset(dat, WTMEC2YR>0), nest=TRUE)
	dsn
}

surveyNumber <- function(seriesName) {
	if(seriesName == '1999-2000') {
		return(1)
	}else if (seriesName == '2001-2002') {
		return(2)
	} else if(seriesName == '2003-2004') {
		return(3)
	} else if(seriesName == '2005-2006') {
		return(4)
	}
	return(NA)
}

numCommonMeasures <- function(dataTable, var1, var2) {
	return(sum(!is.na(dataTable[, var1]) & !is.na(dataTable[,var2])))
}

varDesc <- function(metaTable, variable) {
	lis <- metaTable[metaTable$var == variable, 'var_desc_ewas']
	return(lis[1])
}

hasZero <- function(arr) {
	sum(!is.na(arr) & arr == 0) > 0
}


computeResiduals <- function(dsn, variable, adjustmentVars) {
	formula <- create_reg_formula(adjustmentVars, variable)
	mod <- svyglm(formula, dsn)
	resid <- residuals(mod)
	toReturn <- mod$survey.design$variables
	toReturn[, variable] <- resid
	return(toReturn)
}

adjustVariable <- function(dat, variable, seriesName, adjustmentVariables, metaTable) {
	##
	adjustBy <- setdiff(adjustmentVariables, variable)
	if(length(adjustBy) == 0) {
		return(dat)
	}
	
	refGroup <- binaryRefGroup(metaTable, variable, seriesName)
	binOrd <- isBinaryOrdinal(metaTable, variable, seriesName)
	
	if(refGroup > 0 | binOrd > 0) {
		return(dat)
	}
	
	dsn <- surveyDesign(dat)
	return(computeResiduals(dsn, variable, adjustBy))
}

saveTable <- function(variable,seriesName, tab, outpath='.') {
	outfilename <- sprintf('%s_%s.Rdata', variable, seriesName)
	filename <- file.path(outpath, outfilename)
	save(variable, tab, file=filename)
}

isBinaryVariable <- function(varname, metaTable) {
	varInfo <- metaTable[metaTable$var == varname, ]
	varInfo <- varInfo[1, ]
	
	if(!is.na(varInfo$is_binary) & varInfo$is_binary == 1) {
		return(TRUE)
	}
	return(FALSE)
}

logvariable <- function(varname, metaTable) {
	## log the biochemistry variables?
	## log the environmental
	## log the nutrients
	
	## first check if the categorical or the binary is set
	
	varInfo <- metaTable[metaTable$var == varname, ]
	varInfo <- varInfo[1, ]
	
	if(!is.na(varInfo$is_binary) & varInfo$is_binary == 1) {
		return(FALSE)
	}
	
	if(!is.na(varInfo$categorical_ref_group)) {
		return(FALSE)
	}
	
	if(!is.na(varInfo$is_ordinal) & varInfo$is_ordinal == 1) {
		return(FALSE)
	}
	
	if(length(grep('^LB', varname))) {
		return(TRUE)
	} else if(length(grep('^UR', varname))) {
		return(TRUE)
	} else if(length(grep('^DR1', varname))) {
		return(TRUE)
	}
	
	return(FALSE)
}

continuousVariable <- function(varname, metaTable) {
	### check if continuous
	varInfo <- metaTable[metaTable$var == varname, ]
	varInfo <- varInfo[1, ]
	if(!is.na(varInfo$is_binary) & varInfo$is_binary == 1) {
		return(FALSE)
	}
	if(!is.na(varInfo$categorical_ref_group)) {
		return(FALSE)
	}
	return(TRUE)
}

stratified_tab <- function(tab, strata_name) {
	switch(strata_name, 
		black_male_18 = subset(tab, black == 1 & RIDAGEYR <= 18 & RIAGENDR == 1),
		black_male_19_30 = subset(tab, black == 1 & RIDAGEYR > 18 & RIDAGEYR <= 30 & RIAGENDR == 1),
		black_male_31_60 = subset(tab, black == 1 & RIDAGEYR > 30 & RIDAGEYR <= 60 & RIAGENDR == 1),
		black_male_60 = subset(tab, black == 1 & RIDAGEYR > 60 &  RIAGENDR == 1),
		
		black_female_18 = subset(tab, black == 1 & RIDAGEYR <= 18 & RIAGENDR == 2),
		black_female_19_30 = subset(tab, black == 1 & RIDAGEYR > 18 & RIDAGEYR <= 30 & RIAGENDR == 2),
		black_female_31_60 = subset(tab, black == 1 & RIDAGEYR > 30 & RIDAGEYR <= 60 & RIAGENDR == 2),
		black_female_60 = subset(tab, black == 1 & RIDAGEYR > 60 &  RIAGENDR == 2),
#		
		mexican_male_18 = subset(tab, mexican == 1 & RIDAGEYR <= 18 & RIAGENDR == 1),
		mexican_male_19_30 = subset(tab, mexican == 1 & RIDAGEYR > 18 & RIDAGEYR <= 30 & RIAGENDR == 1),
		mexican_male_31_60 = subset(tab, mexican == 1 & RIDAGEYR > 30 & RIDAGEYR <= 60 & RIAGENDR == 1),
		mexican_male_60 = subset(tab, mexican == 1 & RIDAGEYR > 60 &  RIAGENDR == 1),
		
		mexican_female_18 = subset(tab, mexican == 1 & RIDAGEYR <= 18 & RIAGENDR == 2),
		mexican_female_19_30 = subset(tab, mexican == 1 & RIDAGEYR > 18 & RIDAGEYR <= 30 & RIAGENDR == 2),
		mexican_female_31_60 = subset(tab, mexican == 1 & RIDAGEYR > 30 & RIDAGEYR <= 60 & RIAGENDR == 2),
		mexican_female_60 = subset(tab, mexican == 1 & RIDAGEYR > 60 &  RIAGENDR == 2),
#
		white_male_18 = subset(tab, white == 1 & RIDAGEYR <= 18 & RIAGENDR == 1),
		white_male_19_30 = subset(tab, white == 1 & RIDAGEYR > 18 & RIDAGEYR <= 30 & RIAGENDR == 1),
		white_male_31_60 = subset(tab, white == 1 & RIDAGEYR > 30 & RIDAGEYR <= 60 & RIAGENDR == 1),
		white_male_60 = subset(tab, white == 1 & RIDAGEYR > 60 &  RIAGENDR == 1),
		
		white_female_18 = subset(tab, white == 1 & RIDAGEYR <= 18 & RIAGENDR == 2),
		white_female_19_30 = subset(tab, white == 1 & RIDAGEYR > 18 & RIDAGEYR <= 30 & RIAGENDR == 2),
		white_female_31_60 = subset(tab, white == 1 & RIDAGEYR > 30 & RIDAGEYR <= 60 & RIAGENDR == 2),
		white_female_60 = subset(tab, white == 1 & RIDAGEYR > 60 &  RIAGENDR == 2),		
		)
}

strata <- c(
	'black_male_18',
	'black_male_19_30',
	'black_male_31_60',
	'black_male_60',
	
	'black_female_18',
	'black_female_19_30',
	'black_female_31_60',
	'black_female_60',
#		
	'mexican_male_18',
	'mexican_male_19_30',
	'mexican_male_31_60',
	'mexican_male_60',
		
	'mexican_female_18',
	'mexican_female_19_30',
	'mexican_female_31_60',
	'mexican_female_60',
#
	'white_male_18',
	'white_male_19_30',
	'white_male_31_60',
	'white_male_60',
		
	'white_female_18',
	'white_female_19_30',
	'white_female_31_60',
	'white_female_60'
)
####
SURVEY <- surveyNumber(SERIESNAME)

if(SURVEY == 1) {
	tab <- bt_a
} else if(SURVEY == 2) {
	tab <- bt_b
} else if(SURVEY == 3) {
	tab <- bt_c
} else {
	tab <- bt_d
}

metaTable <- merge(metaTable, data.frame(var=colnames(tab)))
otherVars <- unique(metaTable$var)
#otherVars <- c('LBXSCR', 'LBXBCD', 'MONTELUKAST')

prepareData <- function(dataFrame, variable, metaTable) {
	cat('preparing data\n')
	if(logvariable(variable, metaTable)) {
		dataFrame[, variable] <- log(dataFrame[, variable] + 0.1)
	}
	for(ii in 1:length(otherVars)) {
		if(otherVars[ii] == variable) {
			next
		}
		if(ii %% 10 == 0) {
			cat(sprintf('%i/%i\n', ii, length(otherVars)))
		}
		if(logvariable(otherVars[ii], metaTable)) {
			dataFrame[, otherVars[ii]] <- log(dataFrame[, otherVars[ii]] + 0.1)
		}
		
	}
	cat('preparing data done\n')
	return(dataFrame)
}

runCorrelations <- function (dataFrame, variable, metaTable, permute=F, prepare=T) {
	## run correlations for a data set
	dat <- dataFrame
	if(prepare) {
		dat <- prepareData(dataFrame, variable, metaTable)
	}
	
	if(permute) {
		dat[, variable] <- sample(dat[, variable])
	}
	dsn <- surveyDesign(dat)
	outFrame <- data.frame()
	for(ii in 1:length(otherVars)) {
		if(otherVars[ii] == variable) {
			next
		}
		if(ii %% 10 == 0) {
			cat(sprintf('%i/%i\n', ii, length(otherVars)))
		}
		numcommon <- numCommonMeasures(dataFrame, variable, otherVars[ii])
		if(numcommon <= 10) {
			outFrame <- rbind(outFrame
				, data.frame(var1=variable, var2=otherVars[ii], SDDSRVYR=SURVEY, N=numcommon, corr=NA, pvalue=NA))
			cat(sprintf('%s,NE\n',otherVars[ii]))
			next
		}
		cat(sprintf('%s',otherVars[ii]))
		crr <- tryCatch(svycorr(variable, otherVars[ii], dsn), error = function(e){list(cor=NA, pvalue=NA)})
		#outFrame <- rbind(outFrame, data.frame(var1=variable, var2=otherVars[ii], SDDSRVYR=SURVEY, N=numcommon, corr=crr$cor, pvalue=crr$pvalue))
		cat(sprintf(',%f,%f\n', crr$cor, crr$pvalue))
	}
	return(outFrame)
}


runPermuteCorrelations <- function (dataFrame, variable, metaTable, numpermute) {
	outFrame <- c()
	dataFrame <- dataFrame[!is.na(dataFrame[, variable]), ]
	dataFrame <- prepareData(dataFrame, variable, metaTable)
	for(ii in 1:numpermute) {
		dataFrame[, variable] <- sample(dataFrame[, variable])
		newFrame <- runCorrelations(dataFrame, variable, metaTable, permute=T, prepare=F)
		newFrame$permute_index <- ii
		outFrame <- rbind(outFrame, newFrame)
	}
	return(outFrame)
}


corrData <- NULL
if(PERMUTE==0) {
	if(!is.null(STRATA)) {
		corrData <- data.frame()
		for(ii in 1:length(strata)) {
			cat('strata:%s\n', strata[ii])
			corrS <- runCorrelations(stratified_tab(tab, strata[ii]), variable, metaTable)
			corrS$strata <- strata[ii]
			corrData <- rbind(corrData, corrS)
		}
	} else {
		corrData <- runCorrelations(tab, variable, metaTable)
	}
} else {
	
	if(!is.null(STRATA)) {
		corrData <- data.frame()
		for(ii in 1:length(strata)) {
			cat('strata:%s\n', strata[ii])
			corrS <- runPermuteCorrelations(stratified_tab(tab, strata[ii]), variable, metaTable,PERMUTE)
			corrS$strata <- strata[ii]
			corrData <- rbind(corrData, corrS)
		}
	} else {
		corrData <- runPermuteCorrelations(tab, variable, metaTable, PERMUTE)
	}
}


save(variable, corrData,PERMUTE, PARAMETRIC, SERIESNAME, DATASCRIPT, file=filename)
