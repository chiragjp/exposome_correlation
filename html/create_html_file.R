# Chirag Patel
#### script to create a html file to browse the circos plots
#### 12/29/14


TITLE <- 'Exposome Globe Browser | Patel Group'
fileout <- 'index.html'
RELPATH_TO_IMG <- '../result/circos/img/globes'


#### load in the result data
load('../result/ewas_results/ewas_clinical_phenome.Rdata')
ewasData.quant <- metaEwas
phenoTable <- unique(ewasData.quant[, c('depvar_c', 'clinical_category', 'var_desc_depvar')])
phenoTable <- phenoTable[order(phenoTable$clinical_category, phenoTable$var_desc_depvar), ]
phenoTable <- rbind(phenoTable, data.frame(depvar_c='diabetes', clinical_category='disease-related outcomes', var_desc_depvar='type 2 diabetes (fasting blood glucose >= 125 mg/dL)'))
phenoTable <- rbind(phenoTable, data.frame(depvar_c='mortality', clinical_category='disease-related outcomes', var_desc_depvar='all-cause mortality (time to death)'))
phenoTable <- rbind(phenoTable, data.frame(depvar_c='preterm', clinical_category='disease-related outcomes', var_desc_depvar='Moms with preterm birth'))
####

individual_globe_html <- function(imgpath, name) {
	src <- sprintf('<html><head><title>Exposome Globe for %s</title>', name)
	src <- paste(src, '<link rel=\"stylesheet\" href=\"../main.css\" />')
	src <- paste(src, '<link href=\'http://fonts.googleapis.com/css?family=Open+Sans\' rel=\'stylesheet\' type=\'text/css\'/>')
	src <- paste(src, '</head>')
	src <- paste(src, '<body>')
	src <- paste(src, '<h3>', name, '</h3>')
	img <- sprintf('<img src=%s width=1000px>', imgpath)
	src <- paste(src, img)
	src <- paste(src, '<br><br>')
	description <- sprintf('Association -log10(p-values) from EWAS are shown as a separate track above each exposure (red points denote EWAS-replicated associations with positive association size and blue points indicate an EWAS-replicated negative association size). Replicated EWAS associations for %s are offset in labeled in red or blue text. Only \"first-degree\" correlations (correlations for validated EWAS findings) are displayed in the globes and displayed in black text.', name)
	src <- paste(src, description)
	src <- paste(src, sprintf('<br><a href=\'../%s\'>Go Back</a></p>', fileout))
	src <- paste(src, '</body></html>')
	return(src)
}

write_globe_html <- function(htmlfilepath, imgpath, name) {
	src <- individual_globe_html(imgpath,name)
	cat(src, file=htmlfilepath)
}

google_analytics <- function() {
	return("<script>
	  (function(i,s,o,g,r,a,m){i[\'GoogleAnalyticsObject\']=r;i[r]=i[r]||function(){
	  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
	  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
	  })(window,document,\'script\',\'//www.google-analytics.com/analytics.js\',\'ga\');

	  ga(\'create\', \'UA-57022462-3\', \'auto\');
	  ga(\'send\', \'pageview\');

	</script>");
}

####


cat('<html>', file=fileout)
cat(sprintf('<head><title>%s</title>', TITLE), file=fileout, append=T)

cat('<link rel=\"stylesheet\" href=\"main.css\" />', file=fileout, append=T)
cat('<link href=\'http://fonts.googleapis.com/css?family=Open+Sans\' rel=\'stylesheet\' type=\'text/css\'/>', file=fileout, append=T)
cat(google_analytics(), file=fileout, append=T)
cat('</head>', file=fileout, append=T)
thumbDirectory <- file.path(RELPATH_TO_IMG, 'thumbnails')
pngDirectory <- file.path(RELPATH_TO_IMG, 'png')

cat('<body>',file=fileout, append=T)
cat('<div id=\"content\">',file=fileout, append=T)
cat('<h2>Phenotype-Exposome Correlation Globe Browser</h2>',file=fileout, append=T)

categories <- unique(phenoTable$clinical_category)

for(ii in 1:length(categories)) {
	subTab <- subset(phenoTable, clinical_category == categories[ii])
	subTab <- subTab[order(subTab$var_desc_depvar), ]
	cat(sprintf('<h3>%s</h3>', categories[ii]),file=fileout, append=T)
	for(jj in 1:nrow(subTab)) {
		filenamePrefix <- subTab[jj, 'depvar_c']
		thumbFile <- sprintf('%s_tn.jpg', filenamePrefix) 
	
		if(!file.exists(file.path(thumbDirectory, thumbFile))) {
			next
		}
	
		pngFile <- sprintf('%s.png', filenamePrefix)
		imgpath <- file.path('..', RELPATH_TO_IMG,'png', pngFile)
		htmlfile <- sprintf('%s.html', filenamePrefix)
		htmlpath <- file.path('globes', htmlfile)
	
		write_globe_html(htmlpath, imgpath, subTab[jj, 'var_desc_depvar'])
	
		img <- sprintf('<img src=%s width=92px>', file.path(RELPATH_TO_IMG, 'thumbnails', thumbFile))
		imgurl <- sprintf("<a href='%s'>%s</a>", htmlpath, img)
		figurehtml <- sprintf('<div class=\'thumbnail\'>%s<br>%s</div>', imgurl, subTab[jj, 'var_desc_depvar'])
		cat(figurehtml, file=fileout, append=T)
		if((jj %% 5) == 0) {
			cat('<br class="clearboth">', file=fileout, append=T)
		}
	}
	cat('<br class="clearboth">', file=fileout, append=T)
}

cat('<br><br>', file=fileout, append=T)
cat('<h3>About the globes</h3>', file=fileout, append=T)
cat('The exposure <a href=\"http://bit.ly/exposomeglobes\">correlation globes</a> above display the correlation between pairs of environmental exposures where at least one of the exposures is associated with a clinical phenotype. Globes are arranged in order of category of clinical phenotype (e.g., body measure parameters, cancer diagnostics, etc.).', file=fileout, append=T)
cat('<br><br>', file=fileout, append=T)
cat('Click on them for greater detail.', file=fileout, append=T)
cat('<br><br>Brought to you by <a href=\'http://www.chiragjpgroup.org\'>chiragjpgroup.org</a>; globes created with <a href=\'http://circos.ca\'>Circos</a> software and code on <a href=\'https://github.com/chiragjp/exposome_correlation\'>GitHub</a>.</html>', file=fileout, append=T)
cat('</div>',file=fileout, append=T)
cat('</body></html>', file=fileout, append=T)