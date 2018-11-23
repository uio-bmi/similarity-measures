
##########################################################################################################################################################################################
## Ad hoc R-Script for creating figure2 in the manuscript titled: 'Beware of the Jaccard: the choice of metric is important and non-trivial in genomic co-localisation analysis'##########
##########################################################################################################################################################################################

#Steps for obtaining the mentioned figure:
# 1- Download the datasets available on-line at ...
# 2- Place them in the same folder where you have created the R-project
# 3- Check that the names of the files are the same as in the script or rename them if you like
# 4- Run the script
# 5- The pdf figure will be created inside the folder of the R-project.

#Clean the R-environment:
rm(list=ls())

# Data with underlying correlation of 0.1
tetra.0.1 = read.table('Galaxy106-[Raw_results__tetrachoric_results,_underlying_correlation_0.1].tabular', header=F)
colnames(tetra.0.1)=c('V1', 'trackName', 'similarity', 'overlap', 'coverage', 'elements')
tetra.0.1$trackID = c(rep('Q',5), rep('R',5))
tetra.0.1$treshold = substr(tetra.0.1$trackName, 49,51)

jaccard.0.1 = read.table('Galaxy105-[Raw_results__jaccard_results,_underlying_correlation_0.1].tabular', header=F)
colnames(jaccard.0.1)=c('V1', 'trackName', 'similarity', 'overlap', 'coverage', 'elements')
jaccard.0.1$trackID = c(rep('Q',5), rep('R',5))
jaccard.0.1$treshold = substr(jaccard.0.1$trackName, 49,51)

forbes.0.1 = read.table('Galaxy104-[Raw_results__forbes_results,_underlying_correlation_0.1].tabular', header=F)
colnames(forbes.0.1)=c('V1', 'trackName', 'similarity', 'overlap', 'coverage', 'elements')
forbes.0.1$trackID = c(rep('Q',5), rep('R',5))
forbes.0.1$treshold = substr(forbes.0.1$trackName, 49,51)

tetra.0.1$normalizedSimilarity = (tetra.0.1$similarity)/(max(tetra.0.1$similarity))
tetra.0.1ForPlot = tetra.0.1[order(as.numeric(tetra.0.1$treshold), decreasing = FALSE), ]
tetra.0.1ForPlot = tetra.0.1ForPlot[tetra.0.1ForPlot$trackID=='R', ]

jaccard.0.1$normalizedSimilarity = (jaccard.0.1$similarity)/(max(jaccard.0.1$similarity) )
jaccard.0.1ForPlot = jaccard.0.1[order(as.numeric(jaccard.0.1$treshold), decreasing = FALSE), ]
jaccard.0.1ForPlot = jaccard.0.1ForPlot[jaccard.0.1ForPlot$trackID=='R', ]

forbes.0.1$normalizedSimilarity = (forbes.0.1$similarity)/(max(forbes.0.1$similarity) )
forbes.0.1ForPlot = forbes.0.1[order(as.numeric(forbes.0.1$treshold), decreasing = FALSE), ]
forbes.0.1ForPlot = forbes.0.1ForPlot[forbes.0.1ForPlot$trackID=='R', ]


# Data with underlying correlation of 0.5
tetra.0.5 = read.table('Galaxy101-[Raw_results__Tetrachoric_results,_underlying_correlation_0.5].tabular', header=F)
colnames(tetra.0.5)=c('V1', 'trackName', 'similarity', 'overlap', 'coverage', 'elements')
tetra.0.5$trackID = c(rep('Q',5), rep('R',5))
tetra.0.5$treshold = substr(tetra.0.5$trackName, 49,51)

jaccard.0.5 = read.table('Galaxy102-[Raw_results__jaccard_results,_underlying_correlation_0.5].tabular', header=F)
colnames(jaccard.0.5)=c('V1', 'trackName', 'similarity', 'overlap', 'coverage', 'elements')
jaccard.0.5$trackID = c(rep('Q',4), rep('R',2), 'Q', rep('R',3))
jaccard.0.5$treshold = substr(jaccard.0.5$trackName, 49,51)

forbes.0.5 = read.table('Galaxy103-[Raw_results__forbes_results,_underlying_correlation_0.5].tabular', header=F)
colnames(forbes.0.5)=c('V1', 'trackName', 'similarity', 'overlap', 'coverage', 'elements')
forbes.0.5$trackID = c(rep('Q',4), rep('R',3), 'Q', rep('R',2))
forbes.0.5$treshold = substr(forbes.0.5$trackName, 49,51)

tetra.0.5$normalizedSimilarity = (tetra.0.5$similarity-0.5 )/(max(tetra.0.5$similarity) -min(tetra.0.5$similarity) )
tetra.0.5ForPlot = tetra.0.5[order(as.numeric(tetra.0.5$treshold), decreasing = FALSE), ]
tetra.0.5ForPlot = tetra.0.5ForPlot[tetra.0.5ForPlot$trackID=='R', ]

jaccard.0.5$normalizedSimilarity = (jaccard.0.5$similarity-0.5)/(max(jaccard.0.5$similarity) -min(jaccard.0.5$similarity) )
jaccard.0.5ForPlot = jaccard.0.5[order(as.numeric(jaccard.0.5$treshold), decreasing = FALSE), ]
jaccard.0.5ForPlot = jaccard.0.5ForPlot[jaccard.0.5ForPlot$trackID=='R', ]

forbes.0.5$normalizedSimilarity = (forbes.0.5$similarity-0.5)/(max(forbes.0.5$similarity) -min(forbes.0.5$similarity) )
forbes.0.5ForPlot = forbes.0.5[order(as.numeric(forbes.0.5$treshold), decreasing = FALSE), ]
forbes.0.5ForPlot = forbes.0.5ForPlot[forbes.0.5ForPlot$trackID=='R', ]


# Data with underlying correlation of 0.7
tetra.0.7 = read.table('Galaxy99-[Raw_results__Tetrachoric_results,_underlying_correlation_0.7].tabular', header=F)
colnames(tetra.0.7)=c('V1', 'trackName', 'similarity', 'overlap', 'coverage', 'elements')
tetra.0.7$trackID = c(rep('Q',5), rep('R',5))
tetra.0.7$treshold = substr(tetra.0.7$trackName, 49,51)

jaccard.0.7 = read.table('Galaxy100-[Raw_results__Jaccard_results,_underlying_correlation_0.7].tabular', header=F)
colnames(jaccard.0.7)=c('V1', 'trackName', 'similarity', 'overlap', 'coverage', 'elements')
jaccard.0.7$trackID = c(rep('Q',3), rep('R',2), 'Q', rep('R',2), 'Q', 'R')
jaccard.0.7$treshold = substr(jaccard.0.7$trackName, 49,51)

forbes.0.7 = read.table('Galaxy98-[Raw_results__Forber_results,_underlying_correlation_0.7].tabular', header=F)
colnames(forbes.0.7)=c('V1', 'trackName', 'similarity', 'overlap', 'coverage', 'elements')
forbes.0.7$trackID = c(rep('Q',3), rep('R',2), 'Q', rep('R',2), 'Q', 'R')
forbes.0.7$treshold = substr(forbes.0.7$trackName, 49,51)


tetra.0.7$normalizedSimilarity = (tetra.0.7$similarity)/max(tetra.0.7$similarity) 
tetra.0.7ForPlot = tetra.0.7[order(as.numeric(tetra.0.7$treshold), decreasing = FALSE), ]
tetra.0.7ForPlot = tetra.0.7ForPlot[tetra.0.7ForPlot$trackID=='R', ]

jaccard.0.7$normalizedSimilarity = (jaccard.0.7$similarity)/max(jaccard.0.7$similarity) 
jaccard.0.7ForPlot = jaccard.0.7[order(as.numeric(jaccard.0.7$treshold), decreasing = FALSE), ]
jaccard.0.7ForPlot = jaccard.0.7ForPlot[jaccard.0.7ForPlot$trackID=='R', ]

forbes.0.7$normalizedSimilarity = (forbes.0.7$similarity)/max(forbes.0.7$similarity) 
forbes.0.7ForPlot = forbes.0.7[order(as.numeric(forbes.0.7$treshold), decreasing = FALSE), ]
forbes.0.7ForPlot = forbes.0.7ForPlot[forbes.0.7ForPlot$trackID=='R', ]


#
tetra.0.7ForPlot$letter = c('a','b','c','d','e'); tetra.0.7ForPlot$color = rep('blue',5)
tetra.0.5ForPlot$letter = c('f','g','h','i','l'); tetra.0.5ForPlot$color = rep('gray',5)
tetra.0.1ForPlot$letter = c('m','n','o','p','q'); tetra.0.1ForPlot$color = rep('green',5)


forbes.0.7ForPlot$letter = c('a','b','c','d','e'); forbes.0.7ForPlot$color = rep('blue',5)
forbes.0.5ForPlot$letter = c('f','g','h','i','l'); forbes.0.5ForPlot$color = rep('gray',5)
forbes.0.1ForPlot$letter = c('m','n','o','p','q'); forbes.0.1ForPlot$color = rep('green',5)

jaccard.0.7ForPlot$letter = c('a','b','c','d','e'); jaccard.0.7ForPlot$color = rep('blue',5)
jaccard.0.5ForPlot$letter = c('f','g','h','i','l'); jaccard.0.5ForPlot$color = rep('gray',5)
jaccard.0.1ForPlot$letter = c('m','n','o','p','q'); jaccard.0.1ForPlot$color = rep('green',5)


tetra=rbind(tetra.0.7ForPlot,tetra.0.5ForPlot,tetra.0.1ForPlot)
tetra=tetra[order(tetra$similarity, decreasing = F), ]

forbes=rbind(forbes.0.7ForPlot,forbes.0.5ForPlot,forbes.0.1ForPlot)
forbes=forbes[order(forbes$similarity, decreasing = F), ]

jaccard=rbind(jaccard.0.7ForPlot,jaccard.0.5ForPlot,jaccard.0.1ForPlot)
jaccard=jaccard[order(jaccard$similarity, decreasing = F), ]
  
  
# Creating the .pdf figure  
pdf('figure2.pdf', width = 6, height = 10)
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







