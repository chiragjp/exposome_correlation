library('getopt')
library('data.table')
spec = matrix(c('directory', 'd', 1, 'character',
				'fileout', 'o', 1, 'character'
				),
				nrow=2, byrow=TRUE);
				
opt <- getopt(spec)
files <- dir(opt$directory, pattern='.Rdata')
corrDatas <- vector("list", length(files))
for(ii in 1:length(files)) {
	print(files[ii])
	load(file.path(opt$directory, files[ii]))
	corrDatas[[ii]] <- corrData
}

### now combine them all
corrDataGather <- rbindlist(corrDatas)
save(corrDataGather, PARAMETRIC, PERMUTE, SERIESNAME, file=opt$fileout)