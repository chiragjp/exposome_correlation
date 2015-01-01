## convcerts the LD file into something for andy
## 11/20/2013
fisher <- function(correl) {
	return(.5 * log((1+correl)/(1-correl)))
}

fisher_var <- function(n) {
	return((1/(n-3)))
}

load('../../LD/corrdata6/entire_2003-2004.Rdata')
series <- '2003-2004'
x1 <- corrDataGather[,1]
x2 <- corrDataGather[,2]
labels <- data.frame(var1=corrDataGather[,1], var2=corrDataGather[,2])
N <- as.numeric(corrDataGather[,4])
correl <- as.numeric(corrDataGather[,5])
fishers <- fisher(as.numeric(corrDataGather[,5]))
fisherSe <- fisher_var(as.numeric(corrDataGather[,4]))
weightedFisher <- fishers*(1/fisherSe)

corDat <- matrix(nrow=nrow(corrDataGather), ncol=5)
corDat[,1] <- correl
corDat[,2] <- N
corDat[,3] <- fishers
corDat[,4] <- fisherSe
corDat[,5] <- weightedFisher

colnames(corDat) <- c('corr', 'N', 'fisher', 'fisher_var', 'weighted_fisher')

save(labels, corDat, series, file='../../LD/corrData6/corrs_deconvo_cjp_2003-2004.Rdata')

