
##########################################################################################################################################################################################
## Ad hoc R-Script for creating Figure C2 in the manuscript titled: 'Beware of the Jaccard: the choice of metric is important and non-trivial in genomic co-localisation analysis'##########
##########################################################################################################################################################################################

#Steps for obtaining the mentioned figure:
# 1- Download the datasets available on-line at https://hyperbrowser.uio.no/sim-measure/ after having imported the history
# 2- Place them in the same folder where you have created the R-project
# 3- Check that the names of the files are the same as in the script or rename them if you like
# 4- Run the script
# 5- The pdf figure will be created inside the folder of the R-project.

#Clean the R-environment:
rm(list=ls())

## Input file names (change the file names accordibgly to your file names)

# Files generated with underlying correlation of 0.1
sorensen.0.1 = read.table('Galaxy113-[Raw_results__Sorensen-Dice_results,_underlying_correlation_0.1].tabular', header=F)
pearson.0.1 = read.table('Galaxy114-[Raw_results__Pearson_correlation_results,_underlying_correlation_0.1].tabular', header=F)

# Files generated with underlying correlation of 0.5
sorensen.0.5 = read.table('Galaxy115-[Raw_results__Sorensen-Dice_results,_underlying_correlation_0.5].tabular', header=F)
pearson.0.5 = read.table('Galaxy116-[Raw_results__Pearson_correlation_results,_underlying_correlation_0.5].tabular', header=F)

# Files generated with underlying correlation of 0.7
sorensen.0.7 = read.table('Galaxy117-[Raw_results__Sorensen-Dice_results,_underlying_correlation_0.7].tabular', header=F)
pearson.0.7 = read.table('Galaxy118-[Raw_results__Pearson_correlation_results,_underlying_correlation_0.7].tabular', header=F)

# Add colnames to the files
nameOfcolumns = c('V1','trackName','similarity', 'overlap', 'coverage', 'elements')

colnames(sorensen.0.1) = nameOfcolumns
colnames(pearson.0.1) = nameOfcolumns
colnames(sorensen.0.5) = nameOfcolumns
colnames(pearson.0.5) = nameOfcolumns
colnames(sorensen.0.7) = nameOfcolumns
colnames(pearson.0.7) = nameOfcolumns

## Create a funtion which extract the TrachID and Treshold from the track name
getTrackInfo = function(dataset, position){
  n = dim(dataset)[1]
  info = 0
  for(i in 1:n){
    info[i] = strsplit(as.character(dataset$trackName), '-')[[i]][position]
  }
  if (position == 7){
    info[info==1] = 'Q'
    info[info==2] = 'R'}
  return(info)
}

options(warn = -1) 

## Extract TrackID and Treshold for each file
sorensen.0.1$trackID = getTrackInfo(sorensen.0.1, position=7)
sorensen.0.1$treshold = getTrackInfo(sorensen.0.1, position=9)

pearson.0.1$trackID = getTrackInfo(pearson.0.1, position=7)
pearson.0.1$treshold = getTrackInfo(pearson.0.1, position=9)

sorensen.0.5$trackID = getTrackInfo(sorensen.0.5, position=7)
sorensen.0.5$treshold = getTrackInfo(sorensen.0.5, position=9)

pearson.0.5$trackID = getTrackInfo(pearson.0.5, position=7)
pearson.0.5$treshold = getTrackInfo(pearson.0.5, position=9)

sorensen.0.7$trackID = getTrackInfo(sorensen.0.7, position=7)
sorensen.0.7$treshold = getTrackInfo(sorensen.0.7, position=9)

pearson.0.7$trackID = getTrackInfo(pearson.0.7, position=7)
pearson.0.7$treshold = getTrackInfo(pearson.0.7, position=9)


## Creating a function which prepares the data for plotting and selecting only the reference tracks (trackID=='R')
preparationDataForPlot = function(dataset, corr){
  newDataset = dataset[order(as.numeric(dataset$treshold), decreasing = FALSE), ]
  newDataset = newDataset[newDataset$trackID=='R', ]
  if(corr==0.7){
    newDataset$letter = c('a','b','c','d','e'); newDataset$color = 'blue'
  }
  if(corr==0.5){
    newDataset$letter = c('f','g','h','i','l'); newDataset$color = 'gray'
  }
  if(corr==0.1){
    newDataset$letter = c('m','n','o','p','q'); newDataset$color = 'green'
  }
  return(newDataset)
}

# Preparing the data for plotting
sorensen.0.1ForPlot = preparationDataForPlot(sorensen.0.1, 0.1)
sorensen.0.5ForPlot = preparationDataForPlot(sorensen.0.5, 0.5)
sorensen.0.7ForPlot = preparationDataForPlot(sorensen.0.7, 0.7)

pearson.0.1ForPlot = preparationDataForPlot(pearson.0.1, 0.1)
pearson.0.5ForPlot = preparationDataForPlot(pearson.0.5, 0.5)
pearson.0.7ForPlot = preparationDataForPlot(pearson.0.7, 0.7)

## Combining the data for each similarity measure and ordering them
sorensen = rbind (sorensen.0.7ForPlot, sorensen.0.5ForPlot, sorensen.0.1ForPlot)
sorensen = sorensen[order(sorensen$similarity, decreasing = F), ]

pearson = rbind(pearson.0.7ForPlot, pearson.0.5ForPlot, pearson.0.1ForPlot)
pearson = pearson[order(pearson$similarity, decreasing = F), ]

  
# Creating the .pdf figure 
figureName = 'figureC2.pdf'

pdf(figureName, width = 6, height = 8)
layout(matrix(c(1,1,1,3,1,1,1,3,2,2,2,4,2,2,2,4,5,5,5,5), ncol=4, byrow=T))
par(mar=c(2,2,2,2), oma=c(3,3,2,2))

### Sorensen 
upTetra=max(max(sorensen.0.1ForPlot$similarity), max(sorensen.0.5ForPlot$similarity), max(sorensen.0.7ForPlot$similarity))
plot(sorensen.0.1ForPlot$treshold, sorensen.0.1ForPlot$similarity, type='b', col='green', 
     ylab='similarity', ylim=c(0, 0.2), axes=F,pch=c('m','n','o','p','q'))
lines(sorensen.0.5ForPlot$treshold, sorensen.0.5ForPlot$similarity, type='b', col='gray', 
      ylab='similarity', pch=c('f','g','h','i','l'))
lines(sorensen.0.7ForPlot$treshold, sorensen.0.7ForPlot$similarity, type='b', col='blue', 
      ylab='similarity', pch=c('a','b','c','d','e'))
axis(2, seq(0,0.2,0.05), seq(0,0.2,0.05))
axis(1, seq(2,4,0.5), c('~2270000','~620000', '~134000','~23000','~3100'), line=1)
mtext('Sørensen-Dice index', side=2, line=3, cex = 0.8)
mtext('elements', side=1, line=3.5, cex = 0.8)

#Pearson
upJaccard=max(max(pearson.0.1ForPlot$similarity), max(pearson.0.5ForPlot$similarity), max(pearson.0.7ForPlot$similarity))
plot(pearson.0.1ForPlot$treshold, pearson.0.1ForPlot$similarity, col='green',pch=c('m','n','o','p','q'), type='b',
      ylab='similarity', ylim=c(0, 0.25), axes=F)
points(pearson.0.5ForPlot$treshold, pearson.0.5ForPlot$similarity, type='b', col='gray', 
      ylab='similarity',pch=c('f','g','h','i','l'))
lines(pearson.0.7ForPlot$treshold, pearson.0.7ForPlot$similarity, type='b', col='blue',pch=c('a','b','c','d','e'), 
      ylab='similarity')

mtext('Pearson Correlation', side=2, line=3, cex = 0.8)

axis(1, seq(2,4,0.5), c('~2270000','~620000', '~134000','~23000','~3100'), line=1)
axis(2, seq(0,0.25,0.05), seq(0,0.25,0.05))

# Plotting the bars
plot(c(1,10), c(1,3),xlim=c(1,3),ylim=c(1,10),pch='',axes=F)
points(rep(2,15), seq(1,10,length=15), pch=sorensen$letter, col=sorensen$color)
mtext('Ranked Sørensen-Dice values', side=2, line=-2, cex = 0.7)

plot(c(1,10), c(1,3),xlim=c(1,3),ylim=c(1,10),pch='',axes=F )
points(rep(2,15), seq(1,10,length=15), pch=pearson$letter, col=pearson$color)
mtext('Ranked Pearson values', side=2, line=-2, cex = 0.7)

# Legend
plot(c(1,2,3,4), c(1,2,3,4), pch='', axes=F)
legend(1,3, legend = c('true underlying correlation=0.7', 'true underlying correlation=0.5', 'true underlying correlation=0.1'), col=c('blue', 'gray', 'green'),
       pch=1, lty=1, box.col = 'white')

dev.off()







