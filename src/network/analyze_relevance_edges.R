library(ggplot2)
source('util.R')
source('relevance_edges.R')

################### plot the CDF
toPlot <- corrData2
toPlot$qvalue_0_or_gt_1 <- 'all correlations'
toPlot[toPlot$num_qvalue_05 > 0, 'qvalue_0_or_gt_1'] <- 'q<=0.05 in at least 1 survey'
toPlot[toPlot$num_qvalue_05 > 1, 'qvalue_0_or_gt_1'] <- 'q<=0.05 in at least 2 surveys [replicated]'
p <- ggplot(toPlot, aes(abs(corr), color=factor(qvalue_0_or_gt_1)))
p <- p + stat_ecdf() + ylab('Cumulative fraction') + xlab('|Correlation|') + theme(legend.position='top', legend.title=element_blank())
p  

################### plot the histogram
toPlot <- subset(corrData2, num_qvalue_05 > 1)
p <- ggplot(toPlot, aes(corr))
p <- p +  stat_bin(aes(y=I(100*..count../sum(..count..))), colour = "blue", fill = "white") 
p <- p + geom_vline(xintercept=0) + xlab('Correlation') + ylab('Percent')  
###################

######### STATS

# average abs(correlation) for all
mean(abs(corrData2$corr))
# MEAN: 0.06
median(abs(corrData2$corr))
# MEDIAN: 0.025
quantile(abs(corrData2$corr))
# IQR: .01, .06
quantile(abs(subset(corrData2, num_qvalue_05 >= 1)$corr))
quantile(abs(subset(corrData2, num_qvalue_05 > 1)$corr))
sum(subset(corrData2, num_qvalue_05 > 1)$corr > 0)
repcorrs <- (subset(corrData2, num_qvalue_05 > 1))$corr

### correlation of correlations between series
# all:
cor.test(corrData2$corr_a, corrData2$corr_b, use='pairwise.complete.obs')
## 
cor.test(corrData2$corr_a, corrData2$corr_c, use='pairwise.complete.obs')
## 
cor.test(corrData2$corr_a, corrData2$corr_d, use='pairwise.complete.obs')
## 
cor.test(corrData2$corr_b, corrData2$corr_c, use='pairwise.complete.obs')
## 
cor.test(corrData2$corr_b, corrData2$corr_d, use='pairwise.complete.obs')
## 
cor.test(corrData2$corr_c, corrData2$corr_d, use='pairwise.complete.obs')
## 

conCord <- cor(corrData2[, c('corr_a', 'corr_b', 'corr_c', 'corr_d')], use='pairwise.complete.obs')

# validated:
cor.test(exposomeNetwork$corr_a, exposomeNetwork$corr_b, use='pairwise.complete.obs')
## 
cor.test(exposomeNetwork$corr_a, exposomeNetwork$corr_c, use='pairwise.complete.obs')
## 
cor.test(exposomeNetwork$corr_a, exposomeNetwork$corr_d, use='pairwise.complete.obs')
## 
cor.test(exposomeNetwork$corr_b, exposomeNetwork$corr_c, use='pairwise.complete.obs')
## 
cor.test(exposomeNetwork$corr_b, exposomeNetwork$corr_d, use='pairwise.complete.obs')
## 
cor.test(exposomeNetwork$corr_c, exposomeNetwork$corr_d, use='pairwise.complete.obs')
## 

###



