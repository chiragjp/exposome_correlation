source('util.R')
load('../../data/corrdata6/bigTable_correlation_nhanes.Rdata')
load('../../data/corrdata6/entire_spearman_2003-2004.Rdata')
# load('../../data/corrdata6/entire_permute_2003-2004.Rdata')
# corrPermute <- corrDataGather;
# > quantile(abs(corrPerm2$corr), na.rm=T, probs=0.95)
#        95% 
# 0.07960532 
# > quantile(abs(corrPerm2$corr), na.rm=T, probs=0.99)
#      99% 
# 0.202113 


varDescEwas <- c('acrylamide', 'allergen test', 'bacterial infection', 'cotinine', 'demographics', 'diakyl', 'dioxins', 'furans', 'heavy metals','hydrocarbons',
'latex', 'melamine', 'nutrients', 'pcbs', 'perchlorate', 'pesticides', 'phenols', 'phthalates', 'phytoestrogens', 'polybrominated ethers',
'polyflourochemicals', 'viral infection', 'volatile compounds', 'pharmaceutical', 'smoking behavior', 'physical fitness', 'food component recall', 'hormone')

metaTable <- metaTable[metaTable$var_desc_ewas %in% varDescEwas & metaTable$series == '2003-2004', ] 
metaTable[metaTable$var == 'LBXCOT', 'var_desc_ewas'] <- 'smoking behavior'
metaTable <- metaTable[metaTable$var != 'income', ]
metaTable[metaTable$var == 'INDFMPIR', 'var_desc_ewas_sub'] <- 'sociodemographics'
metaTable[which(metaTable$var_desc_ewas_sub == 'education'), 'var_desc_ewas_sub'] <- 'sociodemographics'
metaTable[which(metaTable$var_desc_ewas_sub == 'occupation'), 'var_desc_ewas_sub'] <- 'sociodemographics'
vars <- metaTable$var
vars <- vars[-grep('LA$', vars)]

threeChars <- substr(vars, 4, 6)
siVariable <- paste('LBD', threeChars, 'SI', sep="")
vars <- vars[!(vars %in% siVariable)]
vars <- vars[-grep('^SXQ', vars)]
vars <- vars[-grep('^RHQ', vars)]
vars <- vars[-grep('^SMD', vars)]
vars <- vars[-grep('^DRD', vars)]
vars <- vars[-grep('^CVD', vars)]
vars <- c(vars, 'CVDVOMAX', 'CVDESVO2')

metaTable <- metaTable[metaTable$var %in% vars, ]
corrData <- corrDataGather[corrDataGather$var1 %in% metaTable$var & corrDataGather$var2 %in% metaTable$var , ]

#### remove redundant ones
key <- vector('list', nrow(corrData))
for(ii in 1:nrow(corrData)) {
	key[ii] <- paste(sort(t(corrData[ii, c('var1', 'var2')])), collapse=';')
}
corrData$key <- unlist(key)
corrData2 <- unique(corrData[, c('SDDSRVYR', 'N', 'corr', 'key')])
corrData2$var1 <- sapply(strsplit(corrData2$key, ';'), function(elem) {elem[1]})
corrData2$var2 <- sapply(strsplit(corrData2$key, ';'), function(elem) {elem[2]})

nodeInfo <- metaTable[, c('var', 'var_desc', 'var_desc_ewas', 'var_desc_ewas_sub')]
nodeInfo[, 'var_desc_ewas_original'] <- nodeInfo$var_desc_ewas
## combine the sub and regular for pesticides and demographics
pest <- which(!is.na(nodeInfo$var_desc_ewas_sub) & nodeInfo$var_desc_ewas == 'pesticides')
nodeInfo[pest, 'var_desc_ewas'] <- paste(nodeInfo$var_desc_ewas[pest], nodeInfo$var_desc_ewas_sub[pest], sep=' ')
demo <- which(!is.na(nodeInfo$var_desc_ewas_sub) & nodeInfo$var_desc_ewas == 'demographics')
nodeInfo[demo, 'var_desc_ewas'] <- paste(nodeInfo$var_desc_ewas[demo], nodeInfo$var_desc_ewas_sub[demo], sep=' ')

corrData2 <- merge(corrData2, nodeInfo, by.x='var1', by.y='var', all.x=T, suffix=c('', '_1'))
corrData2 <- merge(corrData2, nodeInfo, by.x='var2', by.y='var', all.x=T, suffix=c('', '_2'))

testBigMatr <- getMatrix(corrData2)
symmMatr <- symmetricMatrix(testBigMatr)
symmMatr[is.na(symmMatr)] <- 0
e <- eigen(symmMatr, symmetric=T)
ev <- e$val
M <- nrow(symmMatr)
Meff <- effective_m(M, ev)
#Meff = 523

countPairsClass <- function(corrData, category_name='pcbs', sub_category_name=NULL) {
	### counts the number of M and pairs
	pcb <- subset(corrData, var_desc_ewas == category_name & var_desc_ewas_2 == category_name)
	if(!is.null(sub_category_name)) {
		pcb <- subset(pcb, var_desc_ewas_sub == sub_category_name & var_desc_ewas_sub_2 == sub_category_name)
	}
	pcbMatr <- symmetricMatrix(getMatrix(pcb))
	numPairs <- sum(!is.na(pcbMatr[(upper.tri(pcbMatr))]))
	M <- nrow(pcbMatr)
	vars <- unique(c(pcb$var1, pcb$var2))
	numSerum <- length(grep('^LBX', vars)) + length(grep('^LBD', vars))
	numUrine <- length(grep('^URX', vars)) + length(grep('^URD', vars))
	return(list(numPairs=numPairs, M=M, num_urine=numUrine, num_serum=numSerum))
}

intraClassEffectiveM <- function(corrData, category_name='pcbs', sub_category_name=NULL) {
	pcb <- subset(corrData, var_desc_ewas == category_name & var_desc_ewas_2 == category_name)
	if(!is.null(sub_category_name)) {
		pcb <- subset(pcb, var_desc_ewas_sub == sub_category_name & var_desc_ewas_sub_2 == sub_category_name)
	}
	pcbMatr <- symmetricMatrix(getMatrix(pcb))
	pcbMatr[is.na(pcbMatr)] <- 0
	e <- eigen(pcbMatr, symmetric=T)
	ev <- e$val
	M <- nrow(pcbMatr)
	Meff <- effective_m(M, ev)
	return(list(M=M, Meff=Meff))
}



### get the count of the variables
tab1 <- corrData2[, c('var1', 'var_desc_ewas')]
tab2 <- corrData2[, c('var2', 'var_desc_ewas_2')]
colnames(tab2) <- c('var1', 'var_desc_ewas')
varTable <- unique(rbind(tab1, tab2))
countTable <- as.data.frame(table(varTable$var_desc_ewas))
colnames(countTable) <- c('var_desc_ewas', 'count')
###

countPairs <- data.frame()
varDescEwas <- unique(nodeInfo$var_desc_ewas)
for(ii in 1:length(varDescEwas)) {
	subTab <- subset(metaTable, var_desc_ewas == varDescEwas[ii])
	frm1 <- as.data.frame(countPairsClass(corrData2, varDescEwas[ii]))
	if(frm1$numPairs == 0) {
		frm2 <- data.frame(M=1, Meff=1)
	} else {
		frm2 <- as.data.frame(intraClassEffectiveM(corrData2, varDescEwas[ii]))
	}
	frm1$M <- NULL
	frm <- cbind(frm1, frm2)
	frm$var_desc_ewas <- varDescEwas[ii]
	countPairs <- rbind(countPairs, frm)
}

####
countPairs <- merge(countTable, countPairs)
countPairs$Mdiff <- countPairs$M - countPairs$Meff
###


library(ggplot2)
library(plyr)
intraClassCorr <- data.frame()
for(ii in 1:length(varDescEwas)) {
	Vd <- varDescEwas[ii]
	frm <- subset(corrData2, var_desc_ewas == Vd & var_desc_ewas_2 == var_desc_ewas)
	intraClassCorr <- rbind(intraClassCorr, frm)
}


toPlot <- intraClassCorr[, c('corr', 'var_desc_ewas')]
corrAll <- corrData2[, c('corr', 'var_desc_ewas')]
corrAll$var_desc_ewas <- 'all'
toPlot <- rbind(toPlot, corrAll)

medianCorr <- ddply(intraClassCorr, .(var_desc_ewas), summarize, med=median(abs(corr), na.rm=T))
medianCorr[order(medianCorr$med, decreasing=T), 'var_desc_ewas']
toPlot$var_desc_ewas_order <- factor(toPlot$var_desc_ewas, levels = c(medianCorr[order(medianCorr$med, decreasing=T), 'var_desc_ewas'], 'all'))
gp <- ggplot(toPlot, aes(x = var_desc_ewas_order, y=abs(corr)))
gp <- gp + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=.07) + xlab("Category") + ylab("Abs(correlation)")

######## don't need below, putting in a a table
number_ticks <- function(n) {function(limits) pretty(limits, n)}
countPairs$var_desc_ewas_order <- factor(countPairs$var_desc_ewas, levels = c('melamine', 'pesticides pyrethyroid', 'pesticides chloroacetanilide', medianCorr[order(medianCorr$med, decreasing=T), 'var_desc_ewas']))
plotPairs1 <- countPairs[, c('var_desc_ewas_order', 'M')]
plotPairs1$effective <- F
plotPairs2 <- countPairs[, c('var_desc_ewas_order', 'Meff')]
plotPairs2$effective <- T
colnames(plotPairs2) <- c('var_desc_ewas_order', 'M', 'effective')
plotPairs <- rbind(plotPairs1, plotPairs2)

meffp <- ggplot(plotPairs, aes(x = var_desc_ewas_order, y=M, fill=effective))
meffp <- meffp + geom_bar(stat='identity', position='dodge') + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position='none') + xlab("Category") + ylab("Number")
meffp + scale_fill_manual(values=c('black','red')) + scale_y_continuous(breaks=number_ticks(20))
meffp2 <- ggplot(plotPairs, aes(x=effective, y=M, group=var_desc_ewas_order)) + geom_line() + geom_point() + scale_y_continuous(breaks=number_ticks(20))
meffp2
#########