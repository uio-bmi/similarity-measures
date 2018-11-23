
##########################################################################################################################################################################################
## Ad hoc R-Script for creating figures 1 and S1 in the manuscript titled: 'Beware of the Jaccard: the choice of metric is important and non-trivial in genomic co-localisation analysis'##########
##########################################################################################################################################################################################

#Steps for obtaining the mentioned figure:
# 1- Download the datasets available on-line at ...
# 2- Place them in the same folder where you have created the R-project
# 3- Check that the names of the files are the same as in the script or rename them if you like
# 4- Run the script
# 5- The pdf figure will be created inside the folder of the R-project
# 6- Re-run the script by changing the names of the files to the files generated with 50%FDR, to obtain the supplementary figure (S1).

#Clean the R-environment:
rm(list=ls())


# GATA2-Data
forbesCombinedGata2 = read.table('Galaxy58-[Raw_results__Forbes_Similarity,_GATA2_5_FDR_vs_Query_GATA1].tabular', header=F)
jaccardCombinedGata2 = read.table('Galaxy59-[Raw_results__Jaccard_Similarity,_GATA2_5_FDR_vs_Query_GATA1].tabular', header=F)
tetraCombinedGata2 = read.table('Galaxy60-[Raw_results__Tetrachoric_correlation,_GATA2_5_FDR_vs_Query_GATA1].tabular', header=F)

colnames(forbesCombinedGata2)=c('','trackName','similarity', 'overlap', 'coverage', 'elements')
colnames(jaccardCombinedGata2)=c('','trackName','similarity', 'overlap', 'coverage', 'elements')
colnames(tetraCombinedGata2)=c('','trackName','similarity', 'overlap', 'coverage', 'elements')


# TF-Data
forbesCombinedTF = read.table('Galaxy57-[Raw_results__forbes_similarity,_TF_5_FDR_vs_GATA1].tabular', header=F)
jaccardCombinedTF = read.table('Galaxy56-[Raw_results__jaccard_similarity,_TF_5_FDR_vs_GATA1].tabular', header=F)
tetraCombinedTF = read.table('Galaxy58-[Raw_results__tetrachoric_similarity,_TF_5_FDR_vs_GATA1].tabular', header=F)

colnames(forbesCombinedTF)=c('','trackName','similarity', 'overlap', 'coverage', 'elements')
colnames(jaccardCombinedTF)=c('','trackName','similarity', 'overlap', 'coverage', 'elements')
colnames(tetraCombinedTF)=c('','trackName','similarity', 'overlap', 'coverage', 'elements')

# Histone modification - Data
forbesCombinedHistone = read.table('Galaxy48-[Raw_results__forbes_Similarity_Histone_Modification_5_FDR_vs_GATA1].tabular', header=F)
jaccardCombinedHistone = read.table('Galaxy49-[Raw_results__jaccard_Similarity_Histone_Modification_5_FDR_vs_GATA1].tabular', header=F)
tetraCombinedHistone = read.table('Galaxy50-[Raw_results__tetrachoric_Similarity_Histone_Modification_5_FDR_vs_GATA1].tabular', header=F)

colnames(forbesCombinedHistone)=c('','trackName','similarity', 'overlap', 'coverage', 'elements')
colnames(jaccardCombinedHistone)=c('','trackName','similarity', 'overlap', 'coverage', 'elements')
colnames(tetraCombinedHistone)=c('','trackName','similarity', 'overlap', 'coverage', 'elements')

#
forbesCombinedGata2$prob = substr(forbesCombinedGata2$trackName, start=56, stop=61)
forbesCombinedTF$prob = substr(forbesCombinedTF$trackName, start=102, stop=107)
forbesCombinedHistone$prob = substr(forbesCombinedHistone$trackName, start=105, stop=110)

jaccardCombinedGata2$prob = substr(jaccardCombinedGata2$trackName, start=56, stop=61)
jaccardCombinedTF$prob = substr(jaccardCombinedTF$trackName, start=102, stop=107)
jaccardCombinedHistone$prob = substr(jaccardCombinedHistone$trackName, start=105, stop=110)

tetraCombinedGata2$prob = substr(tetraCombinedGata2$trackName, start=56, stop=61)
tetraCombinedTF$prob = substr(tetraCombinedTF$trackName, start=102, stop=107)
tetraCombinedHistone$prob = substr(tetraCombinedHistone$trackName, start=105, stop=110)


# Be-aware that in this part you will get some warnings, but they are not relevant for our purpose since releted to the
# impossibility of calculating summary statics on not numerical variables

# Summary Forbes  
forbesCombinedGata2_mean = aggregate(forbesCombinedGata2, by=list(forbesCombinedGata2$prob), FUN=mean)
forbesCombinedGata2_mean$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
forbesCombinedGata2_sd = aggregate(forbesCombinedGata2, by=list(forbesCombinedGata2$prob), FUN=sd)
forbesCombinedGata2_sd$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)

forbesCombinedTF_mean = aggregate(forbesCombinedTF, by=list(forbesCombinedTF$prob), FUN=mean)
forbesCombinedTF_mean$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
forbesCombinedTF_sd = aggregate(forbesCombinedTF, by=list(forbesCombinedTF$prob), FUN=sd)
forbesCombinedTF_sd$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)

forbesCombinedHistone_mean = aggregate(forbesCombinedHistone, by=list(forbesCombinedHistone$prob), FUN=mean)
forbesCombinedHistone_mean$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
forbesCombinedHistone_sd = aggregate(forbesCombinedHistone, by=list(forbesCombinedHistone$prob), FUN=sd)
forbesCombinedHistone_sd$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
warnings()

#Summary Jaccard
jaccardCombinedGata2_mean = aggregate(jaccardCombinedGata2, by=list(jaccardCombinedGata2$prob), FUN=mean)
jaccardCombinedGata2_mean$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
jaccardCombinedGata2_sd = aggregate(jaccardCombinedGata2, by=list(jaccardCombinedGata2$prob), FUN=sd)
jaccardCombinedGata2_sd$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)

jaccardCombinedTF_mean = aggregate(jaccardCombinedTF, by=list(jaccardCombinedTF$prob), FUN=mean)
jaccardCombinedTF_mean$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
jaccardCombinedTF_sd = aggregate(jaccardCombinedTF, by=list(jaccardCombinedTF$prob), FUN=sd)
jaccardCombinedTF_sd$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)

jaccardCombinedHistone_mean = aggregate(jaccardCombinedHistone, by=list(jaccardCombinedHistone$prob), FUN=mean)
jaccardCombinedHistone_mean$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
jaccardCombinedHistone_sd = aggregate(jaccardCombinedHistone, by=list(jaccardCombinedHistone$prob), FUN=sd)
jaccardCombinedHistone_sd$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)

#Summary Tetrachoric correlation
tetraCombinedGata2_mean = aggregate(tetraCombinedGata2, by=list(tetraCombinedGata2$prob), FUN=mean)
tetraCombinedGata2_mean$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
tetraCombinedGata2_sd = aggregate(tetraCombinedGata2, by=list(tetraCombinedGata2$prob), FUN=sd)
tetraCombinedGata2_sd$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)

tetraCombinedTF_mean = aggregate(tetraCombinedTF, by=list(tetraCombinedTF$prob), FUN=mean)
tetraCombinedTF_mean$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
tetraCombinedTF_sd = aggregate(tetraCombinedTF, by=list(tetraCombinedTF$prob), FUN=sd)
tetraCombinedTF_sd$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)

tetraCombinedHistone_mean = aggregate(tetraCombinedHistone, by=list(tetraCombinedHistone$prob), FUN=mean)
tetraCombinedHistone_mean$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
tetraCombinedHistone_sd = aggregate(tetraCombinedHistone, by=list(tetraCombinedHistone$prob), FUN=sd)
tetraCombinedHistone_sd$prob=c(0.0001, 0.0002,0.001,0.005,0.03,0.15,0.5,1.0, 0.00005)
######


# Creating the .pdf figure
pdf('figure1.pdf', width = 6, height = 10)
#jpeg('figure1.jpg', width = 6, height = 10, units = "in", res=350)
layout(matrix(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4), ncol=3, byrow=T))

par( mar=c(2,2,4,2), oma=c(1,3,2,2))
plot(1:9, forbesCombinedGata2_mean$similarity[order(forbesCombinedGata2_mean$prob[1:9], decreasing=TRUE)], ylim=c(-150,150), type='b', col='blue', 
     ylab='similarity', xlab='', axes=F)

lines(1:9, forbesCombinedGata2_mean$similarity[order(forbesCombinedGata2_mean$prob[1:9], decreasing = TRUE)] + 
         forbesCombinedGata2_sd$similarity[order(forbesCombinedGata2_sd$prob[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)
lines(1:9, forbesCombinedGata2_mean$similarity[order(forbesCombinedGata2_mean$prob[1:9], decreasing = TRUE)] -
        forbesCombinedGata2_sd$similarity[order(forbesCombinedGata2_sd$prob[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)

lines(1:9, forbesCombinedHistone_mean$similarity[order(forbesCombinedHistone_mean$prob[1:9], decreasing=TRUE)], col='gray', type='b')
lines(1:9, forbesCombinedHistone_mean$similarity[order(forbesCombinedHistone_mean$prob[1:9], decreasing = TRUE)] + 
        forbesCombinedHistone_sd$similarity[order(forbesCombinedHistone_sd$prob[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)
lines(1:9, forbesCombinedHistone_mean$similarity[order(forbesCombinedHistone_mean$prob[1:9], decreasing = TRUE)] -
        forbesCombinedHistone_sd$similarity[order(forbesCombinedHistone_sd$prob[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)

lines(1:9, forbesCombinedTF_mean$similarity[order(forbesCombinedTF_mean$prob[1:9], decreasing=TRUE)], col='green', type='b')
lines(1:9, forbesCombinedTF_mean$similarity[order(forbesCombinedTF_mean$prob[1:9], decreasing = TRUE)] + 
        forbesCombinedTF_sd$similarity[order(forbesCombinedTF_sd$prob[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
lines(1:9, forbesCombinedTF_mean$similarity[order(forbesCombinedTF_mean$prob[1:9], decreasing = TRUE)] -
        forbesCombinedTF_sd$similarity[order(forbesCombinedTF_sd$prob[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
axis(2, seq(-150,150,50), seq(-150,150,50))
axis(1, seq(1,9,1), c('100%', '50%','15%','3%','0.5%','0.1%', '0.02%', '0.01%','0.005%'), line=1)
mtext('Forbes coefficient', side=1, line=-14.5, cex=0.7, adj=0.04)

plot(1:9, jaccardCombinedGata2_mean$similarity[order(jaccardCombinedGata2_mean$prob[1:9], decreasing=TRUE)], ylim=c(0,0.04), type='b', col='blue', 
     ylab='similarity', xlab='', axes=F)

lines(1:9, jaccardCombinedGata2_mean$similarity[order(jaccardCombinedGata2_mean$prob[1:9], decreasing = TRUE)] + 
        jaccardCombinedGata2_sd$similarity[order(jaccardCombinedGata2_sd$prob[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)
lines(1:9, jaccardCombinedGata2_mean$similarity[order(jaccardCombinedGata2_mean$prob[1:9], decreasing = TRUE)] -
        jaccardCombinedGata2_sd$similarity[order(jaccardCombinedGata2_sd$prob[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)

lines(1:9, jaccardCombinedHistone_mean$similarity[order(jaccardCombinedHistone_mean$prob[1:9], decreasing=TRUE)], col='gray', type='b')
lines(1:9, jaccardCombinedHistone_mean$similarity[order(jaccardCombinedHistone_mean$prob[1:9], decreasing = TRUE)] + 
        jaccardCombinedHistone_sd$similarity[order(jaccardCombinedHistone_sd$prob[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)
lines(1:9, jaccardCombinedHistone_mean$similarity[order(jaccardCombinedHistone_mean$prob[1:9], decreasing = TRUE)] -
        jaccardCombinedHistone_sd$similarity[order(jaccardCombinedHistone_sd$prob[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)

lines(1:9, jaccardCombinedTF_mean$similarity[order(jaccardCombinedTF_mean$prob[1:9], decreasing=TRUE)], col='green', type='b')
lines(1:9, jaccardCombinedTF_mean$similarity[order(jaccardCombinedTF_mean$prob[1:9], decreasing = TRUE)] + 
        jaccardCombinedTF_sd$similarity[order(jaccardCombinedTF_sd$prob[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
lines(1:9, jaccardCombinedTF_mean$similarity[order(jaccardCombinedTF_mean$prob[1:9], decreasing = TRUE)] -
        jaccardCombinedTF_sd$similarity[order(jaccardCombinedTF_sd$prob[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
axis(2, seq(0,0.04,0.01), seq(0,0.04,0.01))
axis(1, seq(1,9,1), c('100%', '50%','15%','3%','0.5%','0.1%', '0.02%', '0.01%','0.005%'), line=1)
mtext('Jaccard index', side=1, line=-14.5, cex=0.7, adj=0.04)

plot(1:9, tetraCombinedGata2_mean$similarity[order(tetraCombinedGata2_mean$prob[1:9], decreasing=TRUE)], ylim=c(-1,1), type='b', col='blue', 
     ylab='similarity', xlab='', axes=F)

lines(1:9, tetraCombinedGata2_mean$similarity[order(tetraCombinedGata2_mean$prob[1:9], decreasing = TRUE)] + 
        tetraCombinedGata2_sd$similarity[order(tetraCombinedGata2_sd$prob[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)
lines(1:9, tetraCombinedGata2_mean$similarity[order(tetraCombinedGata2_mean$prob[1:9], decreasing = TRUE)] -
        tetraCombinedGata2_sd$similarity[order(tetraCombinedGata2_sd$prob[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)

lines(1:9, tetraCombinedHistone_mean$similarity[order(tetraCombinedHistone_mean$prob[1:9], decreasing=TRUE)], col='gray', type='b')
lines(1:9, tetraCombinedHistone_mean$similarity[order(tetraCombinedHistone_mean$prob[1:9], decreasing = TRUE)] + 
        tetraCombinedHistone_sd$similarity[order(tetraCombinedHistone_sd$prob[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)
lines(1:9, tetraCombinedHistone_mean$similarity[order(tetraCombinedHistone_mean$prob[1:9], decreasing = TRUE)] -
        tetraCombinedHistone_sd$similarity[order(tetraCombinedHistone_sd$prob[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)

lines(1:9, tetraCombinedTF_mean$similarity[order(tetraCombinedTF_mean$prob[1:9], decreasing=TRUE)], col='green', type='b')
lines(1:9, tetraCombinedTF_mean$similarity[order(tetraCombinedTF_mean$prob[1:9], decreasing = TRUE)] + 
        tetraCombinedTF_sd$similarity[order(tetraCombinedTF_sd$prob[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
lines(1:9, tetraCombinedTF_mean$similarity[order(tetraCombinedTF_mean$prob[1:9], decreasing = TRUE)] -
        tetraCombinedTF_sd$similarity[order(tetraCombinedTF_sd$prob[1:9], decreasing = TRUE)], col='green', type='l', lty=2)

axis(2, seq(-1,1,0.5), seq(-1,1,0.5))
axis(1, seq(1,9,1), c('100%', '50%','15%','3%','0.5%','0.1%', '0.02%', '0.01%','0.005%'), line=1)
mtext('Tetrachoric correlation', side=1, line=-14.5, cex=0.7, adj=0.04)
mtext('% of elements kept per track', side=1.5, line=3.5, cex = 0.7)

# Legend
plot(c(1,2,3,4), c(1,2,3,4), pch='', axes=F)
legend(1,3.5, legend = c('GATA2', 'Histone modification', 'Transcription factor'), col=c('blue', 'gray', 'green'),
       pch=1, lty=1, box.col = 'white')

dev.off()



















