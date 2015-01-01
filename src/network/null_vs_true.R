load('../../data/corrdata6/entire_permute_2003-2004.Rdata')
randomData <- corrDataGather;
load('../../data/corrdata6/entire_spearman_2003-2004.Rdata')

library(ggplot2)
dat <- data.frame(corr=c(randomData$corr, corrDataGather$corr), permute=c(rep(1, length(randomData$corr)), rep(0, length(corrDataGather$corr))))
dat<- dat[complete.cases(dat), ]
toPlot <- dat[sample(1:nrow(dat), 100000), ]
p <- ggplot(dat, aes(abs(corr), color=permute))
p + stat_ecdf()