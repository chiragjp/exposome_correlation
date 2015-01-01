# Chirag J Patel
# utility functions to produce circos plots for exposome correlations
# 12/29/14

######### USE CIRCOS TO PLOT.
source('circos_conf_file.R')
writeConfFile <- function(directory='./img', png_file_name='', link_file_name='', plot_file_name='', plot_text_file_name='', fileout="") {
	cat(circos_conf_file(directory, png_file_name, link_file_name, plot_file_name, plot_text_file_name), file=fileout)
	cat("\n",file=fileout, append=T)
}

writeKaryotypes <- function(varDesc, fileout="") {
	## use chr1 .. chr24
	## order this appropriately
	cat("",file=fileout)
	varDescEwas <- unique(varDesc$var_desc_ewas)
	chrColors <- paste('chr', 1:24, sep="")
	for(ii in 1:length(varDescEwas)) {
		N <- nrow(subset(varDesc, var_desc_ewas == varDescEwas[ii]))
		colr <- chrColors[(ii %% length(chrColors))+1]
		## search and replace spaces with _
		txt <- gsub(" ", "_", varDescEwas[ii])
		ln <- paste('chr -', txt , txt, 0, N, colr)
		cat(ln, file=fileout, append=T)
		cat("\n", file=fileout, append=T)
	}
	
	## now draw the bands
	for(ii in 1:length(varDescEwas)) {
		chr <- gsub(" ", "_", varDescEwas[ii])
		vars <- subset(varDesc, var_desc_ewas==varDescEwas[ii])
		colr <- 'white'
		for(jj in 1:nrow(vars)) {
			txt <- gsub(" ", "_", vars$var_desc[jj])
			x <- vars$x_coord[jj]
			ln <- paste('band', chr, txt,txt, x-1, x, colr)
			cat(ln, file=fileout,append=T)
			cat("\n", file=fileout, append=T)
		}	
	}
	
}

writeLinks <- function(pairsData, varDesc, fileout="") {
	# first line is from
	# second line is to.
	# uniq id, chr, coordinate
	# uniq id, chr, coordinate
	cat("",file=fileout)
	for(ii in 1:nrow(pairsData)) {
		uniqId <- paste(pairsData[ii, 'var1'], pairsData[ii, 'var2'], sep='_')
		varname1 <- pairsData[ii, 'var1']
		varname2 <- pairsData[ii, 'var2']
		chr1 <- pairsData[ii, 'var_desc_ewas1']
		chr2 <- pairsData[ii, 'var_desc_ewas2']
		corr <- pairsData[ii, 'corr']
		key1 <- gsub(" ", "_", chr1)
		key2 <- gsub(" ", "_", chr2)
		x1 <- subset(varDesc, var_desc_ewas == chr1 & var == varname1)$x_coord[1]
		x2 <- subset(varDesc, var_desc_ewas == chr2 & var == varname2)$x_coord[1]
		ln1 <- paste(key1, x1-1, x1)
		ln2 <- paste(key2, x2-1, x2)
		optionString <- sprintf(" corr=%f", corr)
		cat(sprintf("%s %s%s\n", ln1, ln2, optionString),file=fileout,append=T)
	}
	
}
writeEwasPlot <- function(ewasData, varDesc, fileout="") {
	# plots the -log10 pvalue
	#chr start end value options	
	cat("",file=fileout)
	for(ii in 1:nrow(ewasData)) {
		varn <- ewasData[ii, 'varname']
		val <- ewasData[ii, 'log10_pvalue_combined']
		effectSize <- ewasData[ii, 'effect_size']
		validated <- ewasData[ii, 'validated']
		varDescSub <- subset(varDesc, var == varn)
		if(nrow(varDescSub)) {
			chr <- varDescSub$var_desc_ewas[1]
			chr <- gsub(" ", "_", chr)
			start <- as.integer(varDescSub$x_coord[1] -1)
			en <- start + 1
			optionString <- sprintf(' validated=%i,effect_size=%f', validated, effectSize)
			cat(sprintf('%s%s\n', paste(chr, start, en, val), optionString),file=fileout,append=T)
		}
	}
}


writeExposureTextLabels <- function(varDesc, fileout="") {
	cat("",file=fileout)
	for(ii in 1:nrow(varDesc)) {
		varn <- varDesc[ii, 'var']
		chr <- varDesc[ii, 'var_desc_ewas']
		chr <- gsub(" ", "_", chr)
		start <- as.integer(varDesc$x_coord[ii] -1)
		val <- gsub(" ", "_", varDesc$var_desc[ii])
		en <- start + 1
		optionString <- sprintf(' validated=%i', varDesc$validated[ii])
		cat(sprintf('%s%s\n', paste(chr, start, en, val), optionString),file=fileout,append=T)
	}
}

uniqueNodeInfo <- function(metaTable) {
	nodeInfo <- data.frame()
	vars <- unique(metaTable$var)
	for(ii in 1:length(vars)) {
		nodeInfo <- rbind(nodeInfo, subset(metaTable, var == vars[ii])[1, ])
	}
	
	varDescEwas <- unique(nodeInfo$var_desc_ewas)
	nodeInfo$x_coord <- NA
	for(ii in 1:length(varDescEwas)) {
		n <- sum(nodeInfo$var_desc_ewas == varDescEwas[ii])
		nodeInfo[nodeInfo$var_desc_ewas == varDescEwas[ii], 'x_coord'] <- 1:n
	}
	
	return(nodeInfo)
}