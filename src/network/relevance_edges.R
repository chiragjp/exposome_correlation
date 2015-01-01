## load in the data
source('util.R')
load('../../data/corrdata6/entire_spearman_merged.Rdata')
corrPerm_a <- addKey((loadCorrFile('../../data/corrdata6/entire_permute_1999-2000.Rdata')))
corrPerm_b <- addKey((loadCorrFile('../../data/corrdata6/entire_permute_2001-2002.Rdata')))
corrPerm_c <- addKey((loadCorrFile('../../data/corrdata6/entire_permute_2003-2004.Rdata')))
corrPerm_d <- addKey((loadCorrFile('../../data/corrdata6/entire_permute_2005-2006.Rdata')))

varDescEwas <- c('acrylamide', 'allergen test', 'bacterial infection', 'cotinine', 'demographics', 'diakyl', 'dioxins', 'furans', 'heavy metals','hydrocarbons',
'latex', 'melamine', 'nutrients', 'pcbs', 'perchlorate', 'pesticides', 'phenols', 'phthalates', 'phytoestrogens', 'polybrominated ethers',
'polyflourochemicals', 'viral infection', 'volatile compounds', 'pharmaceutical', 'smoking behavior', 'physical fitness', 'food component recall', 'hormone')

metaTable <- metaTable[metaTable$var_desc_ewas %in% varDescEwas , ] 
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
vars <- vars[-grep('^VTQ', vars)]
vars <- vars[-grep('^SMD', vars)]
vars <- vars[-grep('^DRD', vars)]
vars <- vars[-grep('^CVD', vars)]
vars <- vars[!(vars %in% c('URXNAL', 'how_long_estrogen_patch', 'how_long_progestin', 'age_started_birth_control', 'how_long_estrogen_progestin', 'age_stopped_birth_control'))]
vars <- c(vars, 'CVDVOMAX', 'CVDESVO2')
vars <- unique(vars)

metaTable <- metaTable[metaTable$var %in% vars, ]
corrData <- corrData[corrData$var1_c %in% metaTable$var & corrData$var2_c %in% metaTable$var , ]
corrPerm_a <- corrPerm_a[corrPerm_a$var1 %in% metaTable$var & corrPerm_a$var2 %in% metaTable$var, ]
corrPerm_b <- corrPerm_b[corrPerm_b$var1 %in% metaTable$var & corrPerm_b$var2 %in% metaTable$var, ]
corrPerm_c <- corrPerm_c[corrPerm_c$var1 %in% metaTable$var & corrPerm_c$var2 %in% metaTable$var, ]
corrPerm_d <- corrPerm_d[corrPerm_d$var1 %in% metaTable$var & corrPerm_d$var2 %in% metaTable$var, ]

corrData$pvalue_c <- pvalue_empirical(corrData$corr_c, corrPerm_c$corr)
corrData$pvalue_a <- pvalue_empirical(corrData$corr_a, corrPerm_a$corr)
corrData$pvalue_b <- pvalue_empirical(corrData$corr_b, corrPerm_b$corr)
corrData$pvalue_d <- pvalue_empirical(corrData$corr_d, corrPerm_d$corr)

corrData$qvalue_a <- p.adjust(corrData$pvalue_a, method='fdr')
corrData$qvalue_b <- p.adjust(corrData$pvalue_b, method='fdr')
corrData$qvalue_c <- p.adjust(corrData$pvalue_c, method='fdr')
corrData$qvalue_d <- p.adjust(corrData$pvalue_d, method='fdr')

corrData$pvalue_train <- apply(corrData[, c('pvalue_a', 'pvalue_b', 'pvalue_d')], 1, pvalue_meta)
corrData$pvalue_combined <- apply(corrData[, c('pvalue_a', 'pvalue_b','pvalue_c','pvalue_d')], 1, pvalue_meta)
corrData$qvalue_train <- p.adjust(corrData$pvalue_train, method='fdr')
corrData$num_qvalue_05 <- apply(corrData[, c('qvalue_a', 'qvalue_b','qvalue_c','qvalue_d')], 1,
function(x) {
	sum(round(x[!is.na(x)],2) <= 0.05)
})

corrData$num_cohorts <- apply(corrData[, c('corr_a', 'corr_b', 'corr_c', 'corr_d')], 1, function(x) sum(!is.na(x)))
corrData$corr_train <- apply(corrData[, c('corr_a', 'corr_b','corr_d', 'N_a', 'N_b', 'N_d')], 1,
function(x) {
	weighted.mean(x[1:3], x[4:6], na.rm=T)
})

corrData$corr <- apply(corrData[, c('corr_a', 'corr_b','corr_c','corr_d', 'N_a', 'N_b', 'N_c', 'N_d')], 1,
function(x) {
	weighted.mean(x[1:4], x[5:8], na.rm=T)
})

corrData$num_train <- apply(corrData[, c('corr_a', 'corr_b', 'corr_d')], 1, function(x) sum(!is.na(x)))
#corrData2 <- subset(corrData, num_train > 0 & !is.na(corr_c))
corrData2 <- subset(corrData, num_cohorts > 1)
ind.a <- which(!is.na(corrData2$corr_a))
ind.b <- which(!is.na(corrData2$corr_b))
ind.c <- which(!is.na(corrData2$corr_c))
ind.d <- which(!is.na(corrData2$corr_d))


length(unique(c(as.character(corrData2$var1_c[ind.a]), as.character(corrData2$var2_c[ind.a]))))
length(unique(c(as.character(corrData2$var1_c[ind.b]), as.character(corrData2$var2_c[ind.b]))))
length(unique(c(as.character(corrData2$var1_c[ind.c]), as.character(corrData2$var2_c[ind.c]))))
length(unique(c(as.character(corrData2$var1_c[ind.d]), as.character(corrData2$var2_c[ind.d]))))
range(corrData2$N_a[ind.a])
range(corrData2$N_b[ind.b])
range(corrData2$N_c[ind.c])
median(corrData2$N_c[ind.c])
range(corrData2$N_d[ind.d])
median(corrData2$N_d[ind.d])


### compute fdr stats
corrData2$validated <- 0
corrData2$validated[corrData2$num_qvalue_05 > 1] <- 1
exposomeNetwork <- subset(corrData2, validated == 1)


