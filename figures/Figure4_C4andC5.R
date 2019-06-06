##########################################################################################################################################################################################
## Ad hoc R-Script for creating figures 4, C4 and C5 in the manuscript titled: 'Beware of the Jaccard: the choice of metric is important and non-trivial in genomic co-localisation analysis'##########
##########################################################################################################################################################################################

#Steps for obtaining the mentioned figure:
# 1- Download the datasets with ranked values
# 2- Place them in the same folder where you have created the R-project
# 3- Check that the name of the file is the same as in this script
# 4- Run the script
# 5- The 3 pdf figures will be created inside the folder of the R-project

#Clean the R-environment:
rm(list=ls())

dataRanksTopOverlap = read.table('ranksWThr.csv', sep='\t', header=T)
head(dataRanksTopOverlap)


thr=nrow(dataRanksTopOverlap)-100
dataRanksTopOverlap$jaccardForbes = (dataRanksTopOverlap$jaccard>thr)*(dataRanksTopOverlap$forbes>thr)
dataRanksTopOverlap$jaccardTetra = (dataRanksTopOverlap$jaccard>thr)*(dataRanksTopOverlap$tetra>thr)
dataRanksTopOverlap$forbesTetra = (dataRanksTopOverlap$forbes>thr)*(dataRanksTopOverlap$tetra>thr)

pdf('figure4.pdf', width = 6, height = 6)
plot(dataRanksTopOverlap$jaccard, dataRanksTopOverlap$forbes, xlab='Jaccard', ylab='Forbes', xlim=c(0, 37000), ylim=c(0,37000))
dev.off()

pdf('figureC4.pdf', width = 6, height = 6)
plot(dataRanksTopOverlap$jaccard, dataRanksTopOverlap$tetra, xlab='Jaccard', ylab='Tetrachoric Correlation', xlim=c(0, 37000), ylim=c(0,37000))
dev.off()

pdf('figureC5.pdf', width = 6, height = 6)
plot(dataRanksTopOverlap$tetra, dataRanksTopOverlap$forbes,  xlab='Forbes', ylab='Tetrachoric Correlation', xlim=c(0, 37000), ylim=c(0,37000))
dev.off()
