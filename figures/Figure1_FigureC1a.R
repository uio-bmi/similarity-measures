##########################################################################################################################################################################################
## Ad hoc R-Script for creating Figure 1 and Figure C1a in the manuscript titled: 'Beware of the Jaccard: the choice of metric is important and non-trivial in genomic co-localisation analysis'##########
##########################################################################################################################################################################################

#Steps for obtaining the mentioned figure:
# 1- Download the datasets available on-line at https://hyperbrowser.uio.no/sim-measure/ after having imported the history
# 2- Place them in the same folder where you have created the R-project
# 3- Check that the names of the files are the same as in the script or rename them if you like
# 4- Run the script
# 5- The pdf figure will be created inside the folder of the R-project.

# 6- For creating figure C1a, re-run the script by changing the input file names and the name of the figure
# in the pdf generator

#Clean the R-environment:
rm(list=ls())

## Input file names (change the file names accordibgly to your file names)
GATA2_input1 = 'Galaxy102-[Raw_results__Forbes_Similarity,_GATA2_5_FDR_vs_Query_GATA1].tabular'
GATA2_input2 = 'Galaxy103-[Raw_results__Jaccard,_GATA2_5_FDR_vs_Query_GATA1].tabular'
GATA2_input3 = 'Galaxy104-[Raw_results__Tetrachoric,_GATA2_5_FDR_vs_Query_GATA1].tabular'

TF_input1 = 'Galaxy57-[Raw_results__forbes_similarity,_TF_5_FDR_vs_GATA1].tabular'
TF_input2 = 'Galaxy56-[Raw_results__jaccard_similarity,_TF_5_FDR_vs_GATA1].tabular'
TF_input3 = 'Galaxy58-[Raw_results__tetrachoric_similarity,_TF_5_FDR_vs_GATA1].tabular'

HistMod_input1 = 'Galaxy48-[Raw_results__forbes_Similarity_Histone_Modification_5_FDR_vs_GATA1].tabular'
HistMod_input2 = 'Galaxy49-[Raw_results__jaccard_Similarity_Histone_Modification_5_FDR_vs_GATA1].tabular'
HistMod_input3 = 'Galaxy50-[Raw_results__tetrachoric_Similarity_Histone_Modification_5_FDR_vs_GATA1].tabular'
  
## Load GATA2 files for the different similarity index
forbesCombinedGata2 = read.table(GATA2_input1, header=F)
jaccardCombinedGata2 = read.table(GATA2_input2, header=F)
tetraCombinedGata2 = read.table(GATA2_input3, header=F)

## Load TF files for the different similarity index
forbesCombinedTF = read.table(TF_input1, header=F)
jaccardCombinedTF = read.table(TF_input2, header=F)
tetraCombinedTF = read.table(TF_input3, header=F)

## Load Histone modification files for the different similarity index
forbesCombinedHistone = read.table(HistMod_input1, header=F)
jaccardCombinedHistone = read.table(HistMod_input2, header=F)
tetraCombinedHistone = read.table(HistMod_input3, header=F)

# Add colnames to the files
nameOfcolumns = c('','trackName','similarity', 'overlap', 'coverage', 'elements')

colnames(forbesCombinedGata2) = nameOfcolumns
colnames(jaccardCombinedGata2) = nameOfcolumns
colnames(tetraCombinedGata2) = nameOfcolumns
colnames(forbesCombinedTF)  = nameOfcolumns
colnames(jaccardCombinedTF) = nameOfcolumns
colnames(tetraCombinedTF) = nameOfcolumns
colnames(forbesCombinedHistone) = nameOfcolumns
colnames(jaccardCombinedHistone) = nameOfcolumns
colnames(tetraCombinedHistone) = nameOfcolumns

## Create a funtion which extract the probability from the track name
getProbability = function(dataset){
  n = dim(dataset)[1]
  prob = 0
  for(i in 1:n){
    prob[i] = strsplit(as.character(dataset$trackName), '-')[[i]][8]
  }
  prob[prob=='5e']='0.00005'
  return(prob)
}

options(warn = -1) 
forbesCombinedGata2$prob = getProbability(forbesCombinedGata2)
forbesCombinedTF$prob = getProbability(forbesCombinedTF)
forbesCombinedHistone$prob = getProbability(forbesCombinedHistone)

jaccardCombinedGata2$prob = getProbability(jaccardCombinedGata2)
jaccardCombinedTF$prob = getProbability(jaccardCombinedTF)
jaccardCombinedHistone$prob = getProbability(jaccardCombinedHistone)

tetraCombinedGata2$prob = getProbability(tetraCombinedGata2)
tetraCombinedTF$prob = getProbability(tetraCombinedTF)
tetraCombinedHistone$prob = getProbability(tetraCombinedHistone)

##Calculate mean and standard deviation of the similarity index for the 10 simulated tracks 
##for each underlying dataset and probability

##Summary forbes
forbesCombinedGata2_mean = aggregate(forbesCombinedGata2, by=list(forbesCombinedGata2$prob), FUN=mean)
forbesCombinedGata2_sd = aggregate(forbesCombinedGata2, by=list(forbesCombinedGata2$prob), FUN=sd)

forbesCombinedTF_mean = aggregate(forbesCombinedTF, by=list(forbesCombinedTF$prob), FUN=mean)
forbesCombinedTF_sd = aggregate(forbesCombinedTF, by=list(forbesCombinedTF$prob), FUN=sd)

forbesCombinedHistone_mean = aggregate(forbesCombinedHistone, by=list(forbesCombinedHistone$prob), FUN=mean)
forbesCombinedHistone_sd = aggregate(forbesCombinedHistone, by=list(forbesCombinedHistone$prob), FUN=sd)

##Summary jaccard
jaccardCombinedGata2_mean = aggregate(jaccardCombinedGata2, by=list(jaccardCombinedGata2$prob), FUN=mean)
jaccardCombinedGata2_sd = aggregate(jaccardCombinedGata2, by=list(jaccardCombinedGata2$prob), FUN=sd)

jaccardCombinedTF_mean = aggregate(jaccardCombinedTF, by=list(jaccardCombinedTF$prob), FUN=mean)
jaccardCombinedTF_sd = aggregate(jaccardCombinedTF, by=list(jaccardCombinedTF$prob), FUN=sd)

jaccardCombinedHistone_mean = aggregate(jaccardCombinedHistone, by=list(jaccardCombinedHistone$prob), FUN=mean)
jaccardCombinedHistone_sd = aggregate(jaccardCombinedHistone, by=list(jaccardCombinedHistone$prob), FUN=sd)

##Summary Tetrachoric correlation
tetraCombinedGata2_mean = aggregate(tetraCombinedGata2, by=list(tetraCombinedGata2$prob), FUN=mean)
tetraCombinedGata2_sd = aggregate(tetraCombinedGata2, by=list(tetraCombinedGata2$prob), FUN=sd)

tetraCombinedTF_mean = aggregate(tetraCombinedTF, by=list(tetraCombinedTF$prob), FUN=mean)
tetraCombinedTF_sd = aggregate(tetraCombinedTF, by=list(tetraCombinedTF$prob), FUN=sd)

tetraCombinedHistone_mean = aggregate(tetraCombinedHistone, by=list(tetraCombinedHistone$prob), FUN=mean)
tetraCombinedHistone_sd = aggregate(tetraCombinedHistone, by=list(tetraCombinedHistone$prob), FUN=sd)


### Plot Figure 1 (remember to change the figure name if you are running the script for Figure C1a)

figureName = 'figure1.pdf'
pdf(figureName, width = 6, height = 10)
layout(matrix(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4), ncol=3, byrow=T))

par( mar=c(2,2,4,2), oma=c(1,3,2,2))
plot(1:9, forbesCombinedGata2_mean$similarity[order(forbesCombinedGata2_mean$Group.1[1:9], decreasing=TRUE)], ylim=c(-50,450), type='b', col='blue', 
     ylab='similarity', xlab='', axes=F)
lines(1:9, forbesCombinedGata2_mean$similarity[order(forbesCombinedGata2_mean$Group.1[1:9], decreasing = TRUE)] + 
                                                       forbesCombinedGata2_sd$similarity[order(forbesCombinedGata2_sd$Group.1[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)
lines(1:9, forbesCombinedGata2_mean$similarity[order(forbesCombinedGata2_mean$Group.1[1:9], decreasing = TRUE)] -
                                                       forbesCombinedGata2_sd$similarity[order(forbesCombinedGata2_sd$Group.1[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)

lines(1:9, forbesCombinedHistone_mean$similarity[order(forbesCombinedHistone_mean$Group.1[1:9], decreasing=TRUE)], col='gray', type='b')

lines(1:9, forbesCombinedHistone_mean$similarity[order(forbesCombinedHistone_mean$Group.1[1:9], decreasing = TRUE)] + 
        forbesCombinedHistone_sd$similarity[order(forbesCombinedHistone_sd$Group.1[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)
lines(1:9, forbesCombinedHistone_mean$similarity[order(forbesCombinedHistone_mean$Group.1[1:9], decreasing = TRUE)] -
        forbesCombinedHistone_sd$similarity[order(forbesCombinedHistone_sd$Group.1[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)


lines(1:9, forbesCombinedTF_mean$similarity[order(forbesCombinedTF_mean$prob[1:9], decreasing=TRUE)], col='green', type='b')

lines(1:9, forbesCombinedTF_mean$similarity[order(forbesCombinedTF_mean$Group.1[1:9], decreasing = TRUE)] + 
        forbesCombinedTF_sd$similarity[order(forbesCombinedTF_sd$Group.1[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
lines(1:9, forbesCombinedTF_mean$similarity[order(forbesCombinedTF_mean$Group.1[1:9], decreasing = TRUE)] -
        forbesCombinedTF_sd$similarity[order(forbesCombinedTF_sd$Group.1[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
axis(2, seq(-50,450,100), seq(-50,450,100))
axis(1, seq(1,9,1), c('100%', '50%','15%','3%','0.5%','0.1%', '0.02%', '0.01%','0.005%'), line=1)
mtext('Forbes coefficient', side=1, line=-14.5, cex=0.7, adj=0.04)


plot(1:9, jaccardCombinedGata2_mean$similarity[order(jaccardCombinedGata2_mean$Group.1[1:9], decreasing=TRUE)], ylim=c(0,0.15), type='b', col='blue', 
     ylab='similarity', xlab='', axes=F)
lines(1:9, jaccardCombinedGata2_mean$similarity[order(jaccardCombinedGata2_mean$Group.1[1:9], decreasing = TRUE)] + 
        jaccardCombinedGata2_sd$similarity[order(jaccardCombinedGata2_sd$Group.1[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)
lines(1:9, jaccardCombinedGata2_mean$similarity[order(jaccardCombinedGata2_mean$Group.1[1:9], decreasing = TRUE)] -
        jaccardCombinedGata2_sd$similarity[order(jaccardCombinedGata2_sd$Group.1[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)

lines(1:9, jaccardCombinedHistone_mean$similarity[order(jaccardCombinedHistone_mean$Group.1[1:9], decreasing=TRUE)], col='gray', type='b')
lines(1:9, jaccardCombinedHistone_mean$similarity[order(jaccardCombinedHistone_mean$Group.1[1:9], decreasing = TRUE)] + 
        jaccardCombinedHistone_sd$similarity[order(jaccardCombinedHistone_sd$Group.1[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)
lines(1:9, jaccardCombinedHistone_mean$similarity[order(jaccardCombinedHistone_mean$Group.1[1:9], decreasing = TRUE)] -
        jaccardCombinedHistone_sd$similarity[order(jaccardCombinedHistone_sd$Group.1[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)


lines(1:9, jaccardCombinedTF_mean$similarity[order(jaccardCombinedTF_mean$Group.1[1:9], decreasing=TRUE)], col='green', type='b')
lines(1:9, jaccardCombinedTF_mean$similarity[order(jaccardCombinedTF_mean$Group.1[1:9], decreasing = TRUE)] + 
        jaccardCombinedTF_sd$similarity[order(jaccardCombinedTF_sd$Group.1[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
lines(1:9, jaccardCombinedTF_mean$similarity[order(jaccardCombinedTF_mean$Group.1[1:9], decreasing = TRUE)] -
        jaccardCombinedTF_sd$similarity[order(jaccardCombinedTF_sd$Group.1[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
axis(2, seq(0,0.15,0.05), seq(0,0.15,0.05))
axis(1, seq(1,9,1), c('100%', '50%','15%','3%','0.5%','0.1%', '0.02%', '0.01%','0.005%'), line=1)
mtext('Jaccard index', side=1, line=-14.5, cex=0.7, adj=0.04)

plot(1:9, tetraCombinedGata2_mean$similarity[order(tetraCombinedGata2_mean$Group.1[1:9], decreasing=TRUE)], ylim=c(-1,1), type='b', col='blue', 
     ylab='similarity', xlab='', axes=F)
lines(1:9, tetraCombinedGata2_mean$similarity[order(tetraCombinedGata2_mean$Group.1[1:9], decreasing = TRUE)] + 
        tetraCombinedGata2_sd$similarity[order(tetraCombinedGata2_sd$Group.1[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)
lines(1:9, tetraCombinedGata2_mean$similarity[order(tetraCombinedGata2_mean$Group.1[1:9], decreasing = TRUE)] -
        tetraCombinedGata2_sd$similarity[order(tetraCombinedGata2_sd$Group.1[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)

lines(1:9, tetraCombinedHistone_mean$similarity[order(tetraCombinedHistone_mean$Group.1[1:9], decreasing=TRUE)], col='gray', type='b')
lines(1:9, tetraCombinedHistone_mean$similarity[order(tetraCombinedHistone_mean$Group.1[1:9], decreasing = TRUE)] + 
        tetraCombinedHistone_sd$similarity[order(tetraCombinedHistone_sd$Group.1[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)
lines(1:9, tetraCombinedHistone_mean$similarity[order(tetraCombinedHistone_mean$Group.1[1:9], decreasing = TRUE)] -
        tetraCombinedHistone_sd$similarity[order(tetraCombinedHistone_sd$Group.1[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)


lines(1:9, tetraCombinedTF_mean$similarity[order(tetraCombinedTF_mean$Group.1[1:9], decreasing=TRUE)], col='green', type='b')
lines(1:9, tetraCombinedTF_mean$similarity[order(tetraCombinedTF_mean$Group.1[1:9], decreasing = TRUE)] + 
        tetraCombinedTF_sd$similarity[order(tetraCombinedTF_sd$Group.1[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
lines(1:9, tetraCombinedTF_mean$similarity[order(tetraCombinedTF_mean$Group.1[1:9], decreasing = TRUE)] -
        tetraCombinedTF_sd$similarity[order(tetraCombinedTF_sd$Group.1[1:9], decreasing = TRUE)], col='green', type='l', lty=2)

axis(2, seq(-1,1,0.5), seq(-1,1,0.5))
axis(1, seq(1,9,1), c('100%', '50%','15%','3%','0.5%','0.1%', '0.02%', '0.01%','0.005%'), line=1)
mtext('Tetrachoric correlation', side=1, line=-14.5, cex=0.7, adj=0.04)
mtext('% of elements kept per track', side=1.5, line=3.5, cex = 0.7)

## Legend
plot(c(1,2,3,4), c(1,2,3,4), pch='', axes=F)
legend(1,3.5, legend = c('GATA2', 'Histone modification', 'Transcription factor'), col=c('blue', 'gray', 'green'),
       pch=1, lty=1, box.col = 'white')

dev.off()



















