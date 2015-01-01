library('getopt')
spec = matrix(c('directory', 'd', 1, 'character',
				'fileout', 'o', 1, 'character'
				),
				nrow=2, byrow=TRUE);
				
opt <- getopt(spec)
files <- dir(opt$directory, pattern='.Rdata')
corrDatas <- vector("list", length(files))
nrows <- c()
for(ii in 1:length(files)) {
	print(files[ii])
	load(file.path(opt$directory, files[ii]))
	corrDatas[[ii]] <- corrData
	nrows <- c(nrows, nrow(corrData))
}

### now combine them all
corrDataGather <- matrix(nrow=sum(nrows), ncol=ncol(corrDatas[[1]]))
startIndex <- 1
for(ii in 1:length(corrDatas)) {
	ncols <- ncol(corrDatas[[ii]])
	corrDataGather[startIndex:(startIndex+nrows[ii]-1), 1:ncols] <- as.matrix(corrDatas[[ii]])
	startIndex <- startIndex+nrows[ii]
}

save(corrDataGather, PARAMETRIC, PERMUTE, SERIESNAME, file=opt$fileout)
