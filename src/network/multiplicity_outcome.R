source('util.R')
load('../../data/corrdata6/bigTable_correlation_nhanes.Rdata')

metaTable <- metaTable[metaTable$series == '2003-2004', ] 

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

bt_c$diabetes <- NA
bt_c[which(bt_c$LBXGLU >= 126), 'diabetes'] <- 1
bt_c[which(bt_c$LBXGLU < 126), 'diabetes'] <- 0


nVars <- subset(metaTable, var_desc_ewas == 'nutrients')$var
nutrCor <- cor(bt_c[, c('diabetes', nVars)], method='spearman', use='pairwise.complete.obs')
median(abs(nutrCor[1, ]))

nVars <- subset(metaTable, var_desc_ewas == 'nutrients')$var
nutrCor <- cor(bt_c[, c('LBDLDL', nVars)], method='spearman', use='pairwise.complete.obs')
median(abs(nutrCor[1, ]))