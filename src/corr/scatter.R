library('getopt')
spec = matrix(c('outpath', 'o', 1, 'character',
				'series', 's', 1, 'character',
				'permute', 'p', 1, 'integer',
				'parametric', 'c', 1, 'integer',
				'script', 'b', 1, 'character', 
				'strata', 'l', 1, 'character'
				),
				nrow=6, byrow=TRUE);
opt <- getopt(spec)

load('~/projects/exposome_corrs/data/corrdata6/bigtable_correlation_nhanes.Rdata')


if(is.null(opt$series)) {series <- '2003-2004';} else {series <- opt$series;}
if(is.null(opt$permute)) {permute <- 0;} else {permute <- opt$permute;}
if(is.null(opt$parametric)) {pearson <- 1;} else {pearson <- opt$parametric;}

#baseoutdir <- file.path('/home/cjp/cluster/bt_corr', series)
baseoutdir <- file.path('~/exposome_corrs/out', series)
if(is.null(opt$outpath)) {outpath <- baseoutdir;} else {outpath <- file.path(baseoutdir, opt$outpath);}
if(is.null(opt$script)) {datascript <- NULL;} else {datascript <- opt$script;} 
if(is.null(opt$strata)) {strata <- NULL;} else {strata <- opt$strata;} 

pathtocommand <- '~/exposome_corrs/src/corrBigTable.R'


### output a file for each
#header <- 'bsub -q normal '
#jobname <- 'bt_corr'

if(series == '1999-2000') {
	tab <- bt_a
} else if(series == '2001-2002') {
	tab <- bt_b
} else if(series == '2003-2004') {
	tab <- bt_c
} else {
	tab <- bt_d
}

variables <- setdiff(colnames(tab), c('SEQN', 'SDDSRVYR', 'SDMVPSU', 'SDMVSTRA'))
for(ii in 1:length(variables)) {
	variable <- variables[ii]
	#loclJob <- paste(jobname, variable, series, sep="_")
	outfile <- file.path(outpath, sprintf('%s_%s.out', variable, series))
	#cmd <- paste(header, sprintf('-J %s -o %s', loclJob, outfile), sep="")
	cmd <- ""
	rscrip <- paste(cmd, sprintf('Rscript %s -v %s -s \'%s\' -p %i -o %s -c %i', pathtocommand, variable, series, permute, outpath, pearson))
	
	if(!is.null(datascript)) {
		rscrip <- paste(rscrip, sprintf('-b %s', datascript))
	}
	if(!is.null(strata)) {
		rscrip <- paste(rscrip, sprintf('-l %s', strata))
	}
	cat(sprintf('%s\n', rscrip))
}