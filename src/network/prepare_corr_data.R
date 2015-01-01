### prepare data for analysis
source('util.R')
load('../../data/corrdata6/bigTable_correlation_nhanes.Rdata')

corrData_c <- addKey(read.csv('../../data/corrdata6/entire_spearman_2003-2004_nodup.csv'))
corrData_d <- addKey(read.csv('../../data/corrdata6/entire_spearman_2005-2006_nodup.csv'))
corrData_b <- addKey(read.csv('../../data/corrdata6/entire_spearman_2001-2002_nodup.csv'))
corrData_a <- addKey(read.csv('../../data/corrdata6/entire_spearman_1999-2000_nodup.csv'))

mergeCorrs <- function(corrData_a, corrData_b, corrData_c, corrData_d) {
	# merge all the series together for analyses
	corrData <- merge(corrData_c, corrData_a, all.x=T, by.x='key', by.y='key', suffixes=c('_c', '_a'))
	corrData$var1_a <- NULL;
	corrData$var2_a <- NULL;
	corrData$SDDSRVYR_a <- NULL;
	corrData$SDDSRVYR_c <- NULL;
	corrData <- merge(corrData, corrData_b, all.x=T, by.x='key', by.y='key')
	corrData$var1 <- NULL;
	corrData$var2 <- NULL;
	corrData$SDDSRVYR <- NULL;
	corrData$N_b <- corrData$N;
	corrData$N <- NULL
	corrData$corr_b <- corrData$corr;
	corrData$corr <- NULL;
	corrData <- merge(corrData, corrData_d, all.x=T, by.x='key', by.y='key')
	corrData$var1 <- NULL;
	corrData$var2 <- NULL;
	corrData$SDDSRVYR <- NULL;
	corrData$N_d <- corrData$N;
	corrData$N <- NULL
	corrData$corr_d <- corrData$corr;
	corrData$corr <- NULL
	return(corrData)
}

corrData <- mergeCorrs(corrData_a, corrData_b, corrData_c, corrData_d)

save(corrData, metaTable, file='../../data/corrdata6/entire_spearman_merged.Rdata')