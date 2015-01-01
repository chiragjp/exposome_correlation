source('util.R')
load('../../data/corrdata6/bigTable_correlation_nhanes.Rdata')
load('../../data/corrdata6/entire_spearman_2003-2004.Rdata')

varDescEwas <- c('acrylamide', 'allergen test', 'bacterial infection', 'cotinine', 'demographics', 'diakyl', 'dioxins', 'furans', 'heavy metals','hydrocarbons',
'latex', 'melamine', 'nutrients', 'pcbs', 'perchlorate', 'pesticides', 'phenols', 'phthalates', 'phytoestrogens', 'polybrominated ethers',
'polyflourochemicals', 'viral infection', 'volatile compounds')
metaTable <- metaTable[metaTable$var_desc_ewas %in% varDescEwas & metaTable$series == '2003-2004', ] 
vars <- metaTable$var
## need to remove LA, SI
vars <- vars[-grep('LA$', vars)]
threeChars <- substr(vars, 4, 6)
siVariable <- paste('LBD', threeChars, 'SI', sep="")
vars <- vars[!(vars %in% siVariable)]
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

corrData2 <- subset(corrData2, abs(corr) > .2)
####
### write this out for cytoscape
corrData2$direction <- sign(corrData2$corr)
lower <- quantile(corrData2[corrData2$corr < 0, 'corr'])
upper <- quantile(corrData2[corrData2$corr > 0, 'corr'])
corrData2$size <- cut(corrData2$corr, c(lower, upper), labels=F)
##write.csv(corrData2, file='~/projects/EWAS perspective/corr_c.csv', row.names=F, quote=F)
#### now need to color 


nodeInfo <- metaTable[, c('var', 'var_desc', 'var_desc_ewas')]
nodeInfo$category <- 'pollutant'
nodeInfo[nodeInfo$var_desc_ewas == 'demographics', 'category'] <- 'demographics'
nodeInfo[nodeInfo$var_desc_ewas %in% c('nutrients', 'phytoestrogens'), 'category'] <- 'nutrients'
nodeInfo[nodeInfo$var_desc_ewas %in% c('viral infection', 'bacterial infection'), 'category'] <- 'infectious agent'
write.table(nodeInfo, file='~/projects/EWAS perspective/node_info.csv', row.names=F, quote=F, sep=";")


