## Chirag J Patel (cp179)
## 12/29/14

#### creates files formatted for circos
source('util.R')
source('circos.R')
source('relevance_edges.R')
## list of EWAS results
source('ewas_results.R')

### now write the .conf and associated circos files
nodeInfo <- uniqueNodeInfo(metaTable)
nodeInfo <- nodeInfo[nodeInfo$var_desc_ewas != 'pharmaceutical', ] # remove drugs 
nodeInfo$series <- NULL;
nodeInfo$var_desc <- gsub(' \\(.+\\)', '', nodeInfo$var_desc)
nodeInfo$var_desc <- gsub('\\(.+\\)', '', nodeInfo$var_desc)
exposomeNetwork$key <- NULL;
exposomeNetwork$var1 <- exposomeNetwork$var1_c
exposomeNetwork$var2 <- exposomeNetwork$var2_c
exposomeNetwork$var1_c <- NULL
exposomeNetwork$var2_c <- NULL
exposomeLinks <- exposomeNetwork[, c('var1', 'var2', 'corr')]
exposomeLinks <- merge(exposomeLinks, nodeInfo, by.x='var1', by.y='var')
exposomeLinks <- merge(exposomeLinks, nodeInfo, by.x='var2', by.y='var', suffix=c('1', '2'))
exposomeLink.sample <- exposomeLinks[exposomeLinks$var1 %in% c('LBXCOT') | exposomeLinks$var2 %in% c('LBXCOT'), ] 



writeKaryotypes(nodeInfo, fileout='../../result/circos/karyotype.exposome.txt')
writeLinks(exposomeLinks, nodeInfo, fileout='../../result/circos/links/links.exposome.txt')
writeLinks(exposomeLink.sample, nodeInfo, fileout='../../result/circos/links/links.exposome.sample.txt')

##### MORTALITY
ewasData.mortality$validated <- ewasData.mortality$validate
ewasData.mortality$log10_pvalue_combined <- -log10(ewasData.mortality$pvalue_combined)
ewasData.mortality$effect_size <- ewasData.mortality$estimate_combined
ewasData.mortality$log10_pvalue_combined[ewasData.mortality$log10_pvalue_combined>=10] <- 10
ewasData.mortality$varname[ewasData.mortality$varname == 'past_smoking'] <- 'current_past_smoking_1'
ewasData.mortality$varname[ewasData.mortality$varname == 'current_smoking'] <- 'current_past_smoking_2'
validatedVarsUp <- subset(ewasData.mortality, validated == 1 & effect_size > 0)$varname
validatedVarsDown <- subset(ewasData.mortality, validated == 1 & effect_size < 0)$varname
validatedVars <- c(validatedVarsUp, validatedVarsDown)
exposomeLinks.mortality <- exposomeLinks[exposomeLinks$var1 %in% validatedVars | exposomeLinks$var2 %in% validatedVars, ]
### need to subset
writeLinks(exposomeLinks.mortality, nodeInfo, fileout='../../result/circos/links/links.exposome.mortality.txt')
writeEwasPlot(ewasData.mortality, nodeInfo, fileout="../../result/circos/plots/plot.mortality.txt")
## writeExposureTextLabel
vars <- unique(c(validatedVars, as.character(exposomeLinks.mortality$var1), as.character(exposomeLinks.mortality$var2)))
mortalityNodes <- nodeInfo[nodeInfo$var %in% vars, ]
mortalityNodes$validated <- 0
mortalityNodes[mortalityNodes$var %in% validatedVarsUp, 'validated'] <- 1
mortalityNodes[mortalityNodes$var %in% validatedVarsDown, 'validated'] <- -1
writeExposureTextLabels(mortalityNodes, fileout="../../result/circos/texts/text.mortality.txt")
##
writeConfFile(directory='./img', png_file_name='mortality', link_file_name='links.exposome.mortality', plot_file_name='plot.mortality', plot_text_file_name='text.mortality', fileout='../../result/circos/mortality.conf')

##### PRETERM BIRTH
ewasData.preterm$effect_size <- ewasData.preterm$estimate
ewasData.preterm$log10_pvalue_combined <- -log10(ewasData.preterm$pvalue)
validatedVarsUp <- as.character(subset(ewasData.preterm, validated == 1 & effect_size > 0)$varname)
validatedVarsDown <- as.character(subset(ewasData.preterm, validated == 1 & effect_size < 0)$varname)
validatedVars <- c(validatedVarsUp, validatedVarsDown)
exposomeLinks.preterm <- exposomeLinks[exposomeLinks$var1 %in% validatedVars | exposomeLinks$var2 %in% validatedVars, ]
### need to subset
writeLinks(exposomeLinks.preterm, nodeInfo, fileout='../../result/circos/links/links.exposome.preterm.txt')
writeEwasPlot(ewasData.preterm, nodeInfo, fileout="../../result/circos/plots/plot.preterm.txt")
## writeExposureTextLabel
vars <- unique(c(validatedVars,as.character(exposomeLinks.preterm$var1), as.character(exposomeLinks.preterm$var2)))
pretermNodes <- nodeInfo[nodeInfo$var %in% vars, ]
pretermNodes$validated <- 0
pretermNodes[pretermNodes$var %in% validatedVarsUp, 'validated'] <- 1
#diabetesNodes[diabetesNodes$var %in% validatedVarsDown, 'validated'] <- -1
writeExposureTextLabels(pretermNodes, fileout="../../result/circos/texts/text.preterm.txt")
writeConfFile(directory='./img', png_file_name='preterm', link_file_name='links.exposome.preterm', plot_file_name='plot.preterm', plot_text_file_name='text.preterm', fileout='../../result/circos/preterm.conf')


#### DIABETES
ewasData.diabetes$log10_pvalue_combined <- -log10(ewasData.diabetes$pvalue_combined)
validatedVarsUp <- as.character(subset(ewasData.diabetes, validated == 1 & effect_size > 0)$varname)
validatedVarsDown <- as.character(subset(ewasData.diabetes, validated == 1 & effect_size < 0)$varname)
validatedVars <- c(validatedVarsUp, validatedVarsDown)
exposomeLinks.diabetes <- exposomeLinks[exposomeLinks$var1 %in% validatedVars | exposomeLinks$var2 %in% validatedVars, ]
### need to subset
writeLinks(exposomeLinks.diabetes, nodeInfo, fileout='../../result/circos/links/links.exposome.diabetes.txt')
writeEwasPlot(ewasData.diabetes, nodeInfo, fileout="../../result/circos/plots/plot.diabetes.txt")
## writeExposureTextLabel
vars <- unique(c(validatedVars,as.character(exposomeLinks.diabetes$var1), as.character(exposomeLinks.diabetes$var2)))
diabetesNodes <- nodeInfo[nodeInfo$var %in% vars, ]
diabetesNodes$validated <- 0
diabetesNodes[diabetesNodes$var %in% validatedVarsUp, 'validated'] <- 1
diabetesNodes[diabetesNodes$var %in% validatedVarsDown, 'validated'] <- -1
writeExposureTextLabels(diabetesNodes, fileout="../../result/circos/texts/text.diabetes.txt")
writeConfFile(directory='./img', png_file_name='diabetes', link_file_name='links.exposome.diabetes', plot_file_name='plot.diabetes', plot_text_file_name='text.diabetes', fileout='../../result/circos/diabetes.conf')
####

#### 

#### QUANTITATIVE phenotypes
### one plot per quant pheno!
ewasData.quant$log10_pvalue_combined <- -log10(ewasData.quant$pvalue_overall)
ewasData.quant$varname_c[ewasData.quant$varname_c == 'cigarette_smoking'] <- 'current_past_smoking_2'
clinVars <- unique(ewasData.quant$depvar_c)
ewasData.quant$effect_size <- ewasData.quant$estimate_overall
ewasData.quant$validated <- ewasData.quant$validate
ewasData.quant$varname <- ewasData.quant$varname_c
for(ii in 1:length(clinVars)) {
	clinVar <- clinVars[ii]
	eData <- subset(ewasData.quant, depvar_c == clinVar)
	validatedVarsUp <- as.character(subset(eData, validated == 1 & effect_size > 0)$varname)
	validatedVarsDown <- as.character(subset(eData, validated == 1 & effect_size < 0)$varname)
	validatedVars <- c(validatedVarsUp, validatedVarsDown)
	if(length(validatedVars) == 0) {
		next
	}
	exposomeLinks.clin <- exposomeLinks[exposomeLinks$var1 %in% validatedVars | exposomeLinks$var2 %in% validatedVars, ]
	writeLinks(exposomeLinks.clin, nodeInfo, fileout=sprintf('../../result/circos/links/links.exposome.%s.txt', clinVar))
	writeEwasPlot(eData, nodeInfo, fileout=sprintf("../../result/circos/plots/plot.%s.txt", clinVar))
	vars <- unique(c(validatedVars,as.character(exposomeLinks.clin$var1), as.character(exposomeLinks.clin$var2)))
	clinNodes <- nodeInfo[nodeInfo$var %in% vars, ]
	clinNodes$validated <- 0
	clinNodes[clinNodes$var %in% validatedVarsUp, 'validated'] <- 1
	clinNodes[clinNodes$var %in% validatedVarsDown, 'validated'] <- -1
	writeExposureTextLabels(clinNodes, fileout=sprintf("../../result/circos/texts/text.%s.txt", clinVar))
	writeConfFile(directory='./img', png_file_name=clinVar, link_file_name=sprintf('links.exposome.%s', clinVar), plot_file_name=sprintf('plot.%s',clinVar), plot_text_file_name=sprintf('text.%s', clinVar), fileout=sprintf('../../result/circos/%s.conf', clinVar))
}


