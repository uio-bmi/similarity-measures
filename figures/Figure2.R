##########################################################################################################################################################################################
## Ad hoc R-Script for creating Figure 2 in the manuscript titled: 'Beware of the Jaccard: the choice of metric is important and non-trivial in genomic co-localisation analysis'##########
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
tetra.0.1 = read.table('Galaxy106-[Raw_results__tetrachoric_results,_underlying_correlation_0.1].tabular', header=F)
jaccard.0.1 = read.table('Galaxy105-[Raw_results__jaccard_results,_underlying_correlation_0.1].tabular', header=F)
forbes.0.1 = read.table('Galaxy104-[Raw_results__forbes_results,_underlying_correlation_0.1].tabular', header=F)

# Files generated with underlying correlation of 0.5
tetra.0.5 = read.table('Galaxy101-[Raw_results__Tetrachoric_results,_underlying_correlation_0.5].tabular', header=F)
jaccard.0.5 = read.table('Galaxy102-[Raw_results__jaccard_results,_underlying_correlation_0.5].tabular', header=F)
forbes.0.5 = read.table('Galaxy103-[Raw_results__forbes_results,_underlying_correlation_0.5].tabular', header=F)

# Files generated with underlying correlation of 0.7
tetra.0.7 = read.table('Galaxy99-[Raw_results__Tetrachoric_results,_underlying_correlation_0.7].tabular', header=F)
jaccard.0.7 = read.table('Galaxy100-[Raw_results__Jaccard_results,_underlying_correlation_0.7].tabular', header=F)
forbes.0.7 = read.table('Galaxy98-[Raw_results__Forber_results,_underlying_correlation_0.7].tabular', header=F)

# Add colnames to the files
nameOfcolumns = c('V1','trackName','similarity', 'overlap', 'coverage', 'elements')

colnames(tetra.0.1) = nameOfcolumns
colnames(jaccard.0.1) = nameOfcolumns
colnames(forbes.0.1) = nameOfcolumns
colnames(tetra.0.5) = nameOfcolumns
colnames(jaccard.0.5) = nameOfcolumns
colnames(forbes.0.5) = nameOfcolumns
colnames(tetra.0.7) = nameOfcolumns
colnames(jaccard.0.7) = nameOfcolumns
colnames(forbes.0.7) = nameOfcolumns

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
tetra.0.1$trackID = getTrackInfo(tetra.0.1, position=7)
tetra.0.1$treshold = getTrackInfo(tetra.0.1, position=9)

jaccard.0.1$trackID = getTrackInfo(jaccard.0.1, position=7)
jaccard.0.1$treshold = getTrackInfo(jaccard.0.1, position=9)

forbes.0.1$trackID = getTrackInfo(forbes.0.1, position = 7)
forbes.0.1$treshold = getTrackInfo(forbes.0.1, position = 9)

tetra.0.5$trackID = getTrackInfo(tetra.0.5, position = 7)
tetra.0.5$treshold = getTrackInfo(tetra.0.5, position = 9)

jaccard.0.5$trackID = getTrackInfo(jaccard.0.5, position = 7)
jaccard.0.5$treshold = getTrackInfo(jaccard.0.5, position = 9)

forbes.0.5$trackID = getTrackInfo(forbes.0.5, position = 7)
forbes.0.5$treshold = getTrackInfo(forbes.0.5, position = 9)

tetra.0.7$trackID = getTrackInfo(tetra.0.7, position = 7)
tetra.0.7$treshold = getTrackInfo(tetra.0.7, position = 9)

jaccard.0.7$trackID = getTrackInfo(jaccard.0.7, position = 7)
jaccard.0.7$treshold = getTrackInfo(jaccard.0.7, position = 9)

forbes.0.7$trackID = getTrackInfo(forbes.0.7, position = 7)
forbes.0.7$treshold = getTrackInfo(forbes.0.7, position = 9)


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
tetra.0.1ForPlot = preparationDataForPlot(tetra.0.1, 0.1)
tetra.0.5ForPlot = preparationDataForPlot(tetra.0.5, 0.5)
tetra.0.7ForPlot = preparationDataForPlot(tetra.0.7, 0.7)

jaccard.0.1ForPlot = preparationDataForPlot(jaccard.0.1, 0.1)
jaccard.0.5ForPlot = preparationDataForPlot(jaccard.0.5, 0.5)
jaccard.0.7ForPlot = preparationDataForPlot(jaccard.0.7, 0.7)

forbes.0.1ForPlot = preparationDataForPlot(forbes.0.1, 0.1)
forbes.0.5ForPlot = preparationDataForPlot(forbes.0.5, 0.5)
forbes.0.7ForPlot = preparationDataForPlot(forbes.0.7, 0.7)


## Combining the data for each similarity measure and ordering them
tetra = rbind (tetra.0.7ForPlot, tetra.0.5ForPlot, tetra.0.1ForPlot)
tetra = tetra[order(tetra$similarity, decreasing = F), ]

forbes = rbind(forbes.0.7ForPlot, forbes.0.5ForPlot, forbes.0.1ForPlot)
forbes = forbes[order(forbes$similarity, decreasing = F), ]

jaccard = rbind(jaccard.0.7ForPlot, jaccard.0.5ForPlot, jaccard.0.1ForPlot)
jaccard = jaccard[order(jaccard$similarity, decreasing = F), ]
  
  
### Plot Figure 2 

figureName = 'figure2.pdf'
pdf(figureName, width = 6, height = 10)
layout(matrix(c(1,1,1,5,1,1,1,5,2,2,2,6,2,2,2,6,3,3,3,7,3,3,3,7,4,4,4,4), ncol=4, byrow=T))
par(mar=c(2,2,2,2), oma=c(3,3,2,2))

###Forbes
upForbes=log(max(max(forbes.0.1ForPlot$similarity), max(forbes.0.5ForPlot$similarity), max(forbes.0.7ForPlot$similarity)))
plot(forbes.0.1ForPlot$treshold, log(forbes.0.1ForPlot$similarity),  col='green',type='b', 
     ylab='similarity', ylim=c(0, upForbes), axes=F,pch=c('m','n','o','p','q'))
lines(forbes.0.5ForPlot$treshold, log(forbes.0.5ForPlot$similarity), type='b', col='gray', 
      ylab='similarity',pch=c('f','g','h','i','l'))
lines(forbes.0.7ForPlot$treshold, log(forbes.0.7ForPlot$similarity), type='b', col='blue', 
      ylab='similarity', pch=c('a','b','c','d','e'))
axis(2, seq(0,round(upForbes),1), seq(0,round(upForbes),1))
axis(1, seq(2,4,0.5), c('~2270000','~620000', '~134000','~23000','~3100'), line=1)
mtext('Forbes coefficient', side=2, line=3, cex = 0.8)

##Jaccard
upJaccard=max(max(jaccard.0.1ForPlot$similarity), max(jaccard.0.5ForPlot$similarity), max(jaccard.0.7ForPlot$similarity))
plot(jaccard.0.1ForPlot$treshold, jaccard.0.1ForPlot$similarity, col='green',pch=c('m','n','o','p','q'), type='b',
      ylab='similarity', ylim=c(0, upJaccard), axes=F)
points(jaccard.0.5ForPlot$treshold, jaccard.0.5ForPlot$similarity, type='b', col='gray', 
      ylab='similarity',pch=c('f','g','h','i','l'))
lines(jaccard.0.7ForPlot$treshold, jaccard.0.7ForPlot$similarity, type='b', col='blue',pch=c('a','b','c','d','e'), 
      ylab='similarity')
mtext('Jaccard index', side=2, line=3, cex = 0.8)
axis(1, seq(2,4,0.5), c('~2270000','~620000', '~134000','~23000','~3100'), line=1)
axis(2, seq(0,upJaccard,0.01), seq(0,upJaccard,0.01))

### Tetrachoric correlation
upTetra=max(max(tetra.0.1ForPlot$similarity), max(tetra.0.5ForPlot$similarity), max(tetra.0.7ForPlot$similarity))
plot(tetra.0.1ForPlot$treshold, tetra.0.1ForPlot$similarity, type='b', col='green', 
     ylab='similarity', ylim=c(0, upTetra), axes=F,pch=c('m','n','o','p','q'))
lines(tetra.0.5ForPlot$treshold, tetra.0.5ForPlot$similarity, type='b', col='gray', 
      ylab='similarity', pch=c('f','g','h','i','l'))
lines(tetra.0.7ForPlot$treshold, tetra.0.7ForPlot$similarity, type='b', col='blue', 
      ylab='similarity', pch=c('a','b','c','d','e'))
axis(2, seq(0,upTetra,0.1), seq(0,upTetra,0.1))
axis(1, seq(2,4,0.5), c('~2270000','~620000', '~134000','~23000','~3100'), line=1)
mtext('Tetrachoric correlation', side=2, line=3, cex = 0.8)
mtext('elements', side=1, line=3.5, cex = 0.8)

# Legend
plot(c(1,2,3,4), c(1,2,3,4), pch='', axes=F)
legend(1,3, legend = c('true underlying correlation=0.7', 'true underlying correlation=0.5', 'true underlying correlation=0.1'), col=c('blue', 'gray', 'green'),
       pch=1, lty=1, box.col = 'white')

# Plotting the bars
plot(c(1,10), c(1,3),xlim=c(1,3),ylim=c(1,10),pch='',axes=F )
points(rep(2,15), seq(1,10,length=15), pch=forbes$letter, col=forbes$color)
mtext('Ranked forbes values', side=2, line=-2, cex = 0.7)

plot(c(1,10), c(1,3),xlim=c(1,3),ylim=c(1,10),pch='',axes=F )
points(rep(2,15), seq(1,10,length=15), pch=jaccard$letter, col=jaccard$color)
mtext('Ranked jaccard values', side=2, line=-2, cex = 0.7)

plot(c(1,10), c(1,3),xlim=c(1,3),ylim=c(1,10),pch='',axes=F)
points(rep(2,15), seq(1,10,length=15), pch=tetra$letter, col=tetra$color)
mtext('Ranked tetrachoric values', side=2, line=-2, cex = 0.7)

dev.off()







