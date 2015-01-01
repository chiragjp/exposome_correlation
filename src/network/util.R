getMatrix <- function(corrData, columnName='corr') {
	vars <- sort(unique(c(corrData$var1, corrData$var2)))
	corrMatrix <- matrix(NA, ncol=length(vars), nrow=length(vars))
	for(ii in 1:nrow(corrData)) {
		# if((ii %% 1000) == 0) {
		# 	print(ii)
		# }
		x <- which(vars==corrData[ii, 'var1'])
		y <- which(vars==corrData[ii, 'var2'])
		corrMatrix[x,y] <- corrData[ii,columnName]
	}
	colnames(corrMatrix) <- vars
	rownames(corrMatrix) <- vars
	return(corrMatrix)
}

symmetricMatrix <- function(aMatrix) {
	### makes a symmetric matrix from one that just has upper.tri elements
	diag(aMatrix) <- 1
	trans <- t(aMatrix)
	aMatrix[lower.tri(aMatrix)] <- trans[lower.tri(trans)]
	return(aMatrix)
} 

effective_m <- function(M, ev) {
	## caculate the effective number of variables in a matrix
	## M total number of variables
	## ev: eigenvalues of matrix
	return(1 + (M-1)*( 1- (var(ev)/M) ))
}
heatmapColors <- function(numColors=16) {
	c1 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=4/6,end=4.0001/6);
	c2 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=1/6,end=1.0001/6);
	c3 <- c(c1,rev(c2)); 
	return(c3)
}
pvalue_empirical <- function(correlations, randomcorrelations) {
	## two sided p-value for permutation
	randomcorrelations <- randomcorrelations[!is.na(randomcorrelations)]
	randDist <- ecdf(abs(randomcorrelations))
	return(1-randDist(abs(correlations)))
}

pvalue_meta <- function(pvals) {
	## meta analysis of p-value
	pvals <- pvals[!is.na(pvals)]
	if(length(pvals) == 0) {
		return(NA)
	}
	chisq2k <- -2*sum(log(pvals))
	degf <- 2*length(pvals)
	return(pchisq(chisq2k, degf, lower.tail=F))
}

pvalue_rho <- function(correlations, samplesizes) {
	toln <- (correlations+1)/(1-correlations)
	fishers <- .5 * log(toln)
	z <- fishers * sqrt((samplesizes-3)/1.06)
	return(pnorm(abs(z), lower.tail=F)*2)
}
uniqueCorrs <- function(corrData) {
	key <- vector('list', nrow(corrData))
	for(ii in 1:nrow(corrData)) {
		if(ii %% 1000 == 0) {
			print(paste(ii, nrow(corrData)))
		}
		key[ii] <- paste(sort(t(corrData[ii, c('var1', 'var2')])), collapse=';')
	}
	corrData$key <- unlist(key)
	corrData2 <- unique(corrData[, c('SDDSRVYR', 'N', 'corr', 'key')])
	corrData2$var1 <- sapply(strsplit(corrData2$key, ';'), function(elem) {elem[1]})
	corrData2$var2 <- sapply(strsplit(corrData2$key, ';'), function(elem) {elem[2]})
	return(corrData2)
}

loadCorrFile <- function(pathToFile='../../data/corrdata6/entire_spearman_2003-2004.Rdata') {
	load(pathToFile)
	return(corrDataGather)
}

cleanUpCorr <- function(corrData) {
	corrData <- as.data.frame(corrData, stringsAsFactors = F)
	colnames(corrData) <- c('var1', 'var2', 'SDDSRVYR', 'N', 'corr')
	corrData$var1 <- as.character(corrData$var1)
	corrData$var2 <- as.character(corrData$var2)
	corrData$N <- as.numeric(corrData$N)
	corrData$corr <- as.numeric(corrData$corr)
	return(corrData)
}

addKey <- function(corrData) {
	corrData$key <- paste(corrData$var1, corrData$var2, sep=";")
	return(corrData)
}
