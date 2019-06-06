##########################################################################################################################################################################################
## Ad hoc R-Script for creating Figure C1b and Figure C1c in the manuscript titled: 'Beware of the Jaccard: the choice of metric is important and non-trivial in genomic co-localisation analysis'##########
##########################################################################################################################################################################################

#Steps for obtaining the mentioned figure:
# 1- Download the datasets available on-line at https://hyperbrowser.uio.no/sim-measure/ after having imported the history
# 2- Place them in the same folder where you have created the R-project
# 3- Check that the names of the files are the same as in the script or rename them if you like
# 4- Run the script
# 5- The pdf figure will be created inside the folder of the R-project.

# 6- For creating figure C1c, re-run the script by changing the input file names and the name of the figure
# in the pdf generator

#Clean the R-environment:
rm(list=ls())

## Input file names (change the file names accordibgly to your file names)
GATA2_input1 = 'Galaxy105-[Raw_results__Sorensen-Dice,_GATA2_5_FDR_vs_Query_GATA1].tabular'
GATA2_input2 = 'Galaxy106-[Raw_results__Pearson_correlation,_GATA2_5_FDR_vs_Query_GATA1].tabular'

TF_input1 = 'Galaxy61-[Raw_results__Sorensen-Dice_measure,_TF_5_FDR_vs_GATA1].tabular'
TF_input2 = 'Galaxy62-[Raw_results__Pearson_correlation,_TF_5_FDR_vs_GATA1].tabular'

HistMod_input1 = 'Galaxy53-[Raw_results__Sorensen-Dice_measure,_Histone_Modification_5_FDR_vs_GATA1].tabular'
HistMod_input2 = 'Galaxy54-[Raw_results__Pearson_correlation,_Histone_Modification_5_FDR_vs_GATA1].tabular'


## Load GATA2 files for the different similarity index
sorensenCombinedGata2 = read.table(GATA2_input1, header=F)
pearsonCombinedGata2 = read.table(GATA2_input2, header=F)

## Load TF files for the different similarity index
sorensenCombinedTF = read.table(TF_input1, header=F)
pearsonCombinedTF = read.table(TF_input2, header=F)

## Load Histone modification files for the different similarity index
sorensenCombinedHistone = read.table(HistMod_input1, header=F)
pearsonCombinedHistone = read.table(HistMod_input2, header=F)

# Add colnames to the files
nameOfcolumns = c('','trackName','similarity', 'overlap', 'coverage', 'elements')

colnames(sorensenCombinedGata2) = nameOfcolumns
colnames(pearsonCombinedGata2) = nameOfcolumns
colnames(sorensenCombinedTF) = nameOfcolumns
colnames(pearsonCombinedTF) = nameOfcolumns
colnames(pearsonCombinedHistone) = nameOfcolumns
colnames(sorensenCombinedHistone) = nameOfcolumns

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
sorensenCombinedGata2$prob = getProbability(sorensenCombinedGata2)
sorensenCombinedTF$prob = getProbability(sorensenCombinedTF)
sorensenCombinedHistone$prob = getProbability(sorensenCombinedHistone)

pearsonCombinedGata2$prob = getProbability(pearsonCombinedGata2)
pearsonCombinedTF$prob = getProbability(pearsonCombinedTF)
pearsonCombinedHistone$prob = getProbability(pearsonCombinedHistone)

## Calculate mean and standard deviation of the similarity index for the 10 simulated tracks for 
## each underlying dataset and probability
sorensenCombinedGata2_mean = aggregate(sorensenCombinedGata2, by=list(sorensenCombinedGata2$prob), FUN=mean)
sorensenCombinedGata2_sd = aggregate(sorensenCombinedGata2, by=list(sorensenCombinedGata2$prob), FUN=sd)

sorensenCombinedTF_mean = aggregate(sorensenCombinedTF, by=list(sorensenCombinedTF$prob), FUN=mean)
sorensenCombinedTF_sd = aggregate(sorensenCombinedTF, by=list(sorensenCombinedTF$prob), FUN=sd)

sorensenCombinedHistone_mean = aggregate(sorensenCombinedHistone, by=list(sorensenCombinedHistone$prob), FUN=mean)
sorensenCombinedHistone_sd = aggregate(sorensenCombinedHistone, by=list(sorensenCombinedHistone$prob), FUN=sd)

##
pearsonCombinedGata2_mean = aggregate(pearsonCombinedGata2, by=list(pearsonCombinedGata2$prob), FUN=mean)
pearsonCombinedGata2_sd = aggregate(pearsonCombinedGata2, by=list(pearsonCombinedGata2$prob), FUN=sd)

pearsonCombinedTF_mean = aggregate(pearsonCombinedTF, by=list(pearsonCombinedTF$prob), FUN=mean)
pearsonCombinedTF_sd = aggregate(pearsonCombinedTF, by=list(pearsonCombinedTF$prob), FUN=sd)

pearsonCombinedHistone_mean = aggregate(pearsonCombinedHistone, by=list(pearsonCombinedHistone$prob), FUN=mean)
pearsonCombinedHistone_sd = aggregate(pearsonCombinedHistone, by=list(pearsonCombinedHistone$prob), FUN=sd)



### Plot Figure C1b (remember to change the figure name if you are running the script for Figure C1c)

figureName = 'figureC1b.pdf'
pdf(figureName, width = 6, height = 8)
layout(matrix(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3), ncol=3, byrow=T))

par( mar=c(2,2,4,2), oma=c(1,3,2,2))
plot(1:9, sorensenCombinedGata2_mean$similarity[order(sorensenCombinedGata2_mean$Group.1, decreasing=TRUE)], type='b', col='blue', 
     ylab='similarity', xlab='', axes=F, ylim=c(0,0.3))

lines(1:9, sorensenCombinedGata2_mean$similarity[order(sorensenCombinedGata2_mean$Group.1[1:9], decreasing = TRUE)] + 
                                                       sorensenCombinedGata2_sd$similarity[order(sorensenCombinedGata2_sd$Group.1[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)
lines(1:9, sorensenCombinedGata2_mean$similarity[order(sorensenCombinedGata2_mean$Group.1[1:9], decreasing = TRUE)] -
        sorensenCombinedGata2_sd$similarity[order(sorensenCombinedGata2_sd$Group.1[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)

lines(1:9, sorensenCombinedHistone_mean$similarity[order(sorensenCombinedHistone_mean$Group.1[1:9], decreasing=TRUE)], col='gray', type='b')

lines(1:9, sorensenCombinedHistone_mean$similarity[order(sorensenCombinedHistone_mean$Group.1[1:9], decreasing = TRUE)] + 
        sorensenCombinedHistone_sd$similarity[order(sorensenCombinedHistone_sd$Group.1[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)
lines(1:9, sorensenCombinedHistone_mean$similarity[order(sorensenCombinedHistone_mean$Group.1[1:9], decreasing = TRUE)] -
        sorensenCombinedHistone_sd$similarity[order(sorensenCombinedHistone_sd$Group.1[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)


lines(1:9, sorensenCombinedTF_mean$similarity[order(sorensenCombinedTF_mean$Group.1[1:9], decreasing=TRUE)], col='green', type='b')

lines(1:9, sorensenCombinedTF_mean$similarity[order(sorensenCombinedTF_mean$Group.1[1:9], decreasing = TRUE)] + 
        sorensenCombinedTF_sd$similarity[order(sorensenCombinedTF_sd$Group.1[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
lines(1:9, sorensenCombinedTF_mean$similarity[order(sorensenCombinedTF_mean$Group.1[1:9], decreasing = TRUE)] -
        sorensenCombinedTF_sd$similarity[order(sorensenCombinedTF_sd$Group.1[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
axis(2, seq(0,0.3,0.05), seq(0,0.3,0.05))
axis(1, seq(1,9,1), c('100%', '50%','15%','3%','0.5%','0.1%', '0.02%', '0.01%','0.005%'), line=1)
mtext('Sørensen-Dice coefficient', side=1, line=-16.5, cex=0.7, adj=0.1)


plot(1:9, pearsonCombinedGata2_mean$similarity[order(pearsonCombinedGata2_mean$Group.1[1:9], decreasing=TRUE)], ylim=c(-0.01,0.3), type='b', col='blue', 
     ylab='similarity', xlab='', axes=F)

lines(1:9, pearsonCombinedGata2_mean$similarity[order(pearsonCombinedGata2_mean$Group.1[1:9], decreasing = TRUE)] + 
        pearsonCombinedGata2_sd$similarity[order(pearsonCombinedGata2_sd$Group.1[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)
lines(1:9, pearsonCombinedGata2_mean$similarity[order(pearsonCombinedGata2_mean$Group.1[1:9], decreasing = TRUE)] -
        pearsonCombinedGata2_sd$similarity[order(pearsonCombinedGata2_sd$Group.1[1:9], decreasing = TRUE)], col='blue', type='l', lty=2)

lines(1:9, pearsonCombinedHistone_mean$similarity[order(pearsonCombinedHistone_mean$Group.1[1:9], decreasing=TRUE)], col='gray', type='b')
lines(1:9, pearsonCombinedHistone_mean$similarity[order(pearsonCombinedHistone_mean$Group.1[1:9], decreasing = TRUE)] + 
        pearsonCombinedHistone_sd$similarity[order(pearsonCombinedHistone_sd$Group.1[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)
lines(1:9, pearsonCombinedHistone_mean$similarity[order(pearsonCombinedHistone_mean$Group.1[1:9], decreasing = TRUE)] -
        pearsonCombinedHistone_sd$similarity[order(pearsonCombinedHistone_sd$Group.1[1:9], decreasing = TRUE)], col='gray', type='l', lty=2)


lines(1:9, pearsonCombinedTF_mean$similarity[order(pearsonCombinedTF_mean$Group.1[1:9], decreasing=TRUE)], col='green', type='b')
lines(1:9, pearsonCombinedTF_mean$similarity[order(pearsonCombinedTF_mean$Group.1[1:9], decreasing = TRUE)] + 
        pearsonCombinedTF_sd$similarity[order(pearsonCombinedTF_sd$Group.1[1:9], decreasing = TRUE)], col='green', type='l', lty=2)
lines(1:9, pearsonCombinedTF_mean$similarity[order(pearsonCombinedTF_mean$Group.1[1:9], decreasing = TRUE)] -
        pearsonCombinedTF_sd$similarity[order(pearsonCombinedTF_sd$Group.1[1:9], decreasing = TRUE)], col='green', type='l', lty=2)

axis(2, seq(0,0.3,0.05), seq(0,0.3,0.05))
axis(1, seq(1,9,1), c('100%', '50%','15%','3%','0.5%','0.1%', '0.02%', '0.01%','0.005%'), line=1)
mtext('Pearson correlation', side=1, line=-16.5, cex=0.7, adj=0.1)
mtext('% of elements kept per track', side=1.5, line=3.5, cex = 0.7)

## Legend
plot(c(1,2,3,4), c(1,2,3,4), pch='', axes=F)
legend(1,3.5, legend = c('GATA2', 'Histone modification', 'Transcription factor'), col=c('blue', 'gray', 'green'),
       pch=1, lty=1, box.col = 'white')

dev.off()



















