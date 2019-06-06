##########################################################################################################################################################################################
## Ad hoc R-Script for creating Figure 5 in the manuscript titled: 'Beware of the Jaccard: the choice of metric is important and non-trivial in genomic co-localisation analysis'##########
##########################################################################################################################################################################################

#Steps for obtaining the mentioned figure:
# 1- Download the datasets available on-line at https://hyperbrowser.uio.no/sim-measure/ after having imported the specific history
# 2- Place them in the same folder where you have created the R-project
# 3- Check that the names of the files are the same as in the script or rename them if you like
# 4- Run the script
# 5- The pdf figure will be created inside the folder of the R-project.

#Clean the R-environment:
rm(list=ls())

## Input file names (change the file names accordibgly to your file names)
# Data forbes
forbes_input0 = 'Galaxy65-[Raw_results__similarity_forbes_all].tabular'
forbes_input1 = 'Galaxy67-[Raw_results__similarity_forbes-input1].tabular'
forbes_input2 = 'Galaxy68-[Raw_results__similarity_forbes-input2].tabular'
forbes_input3 = 'Galaxy69-[Raw_results__similarity_forbes-input3].tabular'
forbes_input4 = 'Galaxy70-[Raw_results__similarity_forbes-input4].tabular'
forbes_input5 = 'Galaxy71-[Raw_results__similarity_forbes-input5].tabular'
forbes_input6 = 'Galaxy72-[Raw_results__similarity_forbes-input6].tabular'
forbes_input7 = 'Galaxy73-[Raw_results__similarity_forbes-input7].tabular'
forbes_input8 = 'Galaxy74-[Raw_results__similarity_forbes-input8].tabular'
forbes_input9 = 'Galaxy75-[Raw_results__similarity_forbes-input9].tabular'
# Data Jaccard
jaccard_input0 = 'Galaxy66-[Raw_results__similarity_jaccard_all].tabular'
jaccard_input1 = 'Galaxy76-[Raw_results__similarity_jaccard-input1].tabular'
jaccard_input2 = 'Galaxy77-[Raw_results__similarity__jaccard-input2].tabular'
jaccard_input3 = 'Galaxy78-[Raw_results__similarity__jaccard-input3].tabular'
jaccard_input4 = 'Galaxy79-[Raw_results__similarity__jaccard-input4].tabular'
jaccard_input5 = 'Galaxy80-[Raw_results__similarity__jaccard-input5].tabular'
jaccard_input6 = 'Galaxy81-[Raw_results__similarity__jaccard-input6].tabular'
jaccard_input7 = 'Galaxy82-[Raw_results__similarity__jaccard-input7].tabular'
jaccard_input8 = 'Galaxy83-[Raw_results__similarity__jaccard-input8].tabular'
jaccard_input9 = 'Galaxy84-[Raw_results__similarity_to__jaccard-input9].tabular'

## Preparing the data
dataForbesAll = read.table(forbes_input0, head=F, sep='\t')[,c(2,3)]
colnames(dataForbesAll) = c('trackName', 'simAll')

dataForbes = read.table(forbes_input1, head=F, sep='\t')[,c(2,3)]
colnames(dataForbes) = c('trackName', 'sim1')
dataForbesAll = merge(dataForbesAll, dataForbes, by.x='trackName', by.y='trackName')

dataForbes = read.table(forbes_input2, head=F, sep='\t')[,c(2,3)]
colnames(dataForbes) = c('trackName', 'sim2')
dataForbesAll = merge(dataForbesAll, dataForbes, by.x='trackName', by.y='trackName')

dataForbes = read.table(forbes_input3, head=F, sep='\t')[,c(2,3)]
colnames(dataForbes) = c('trackName', 'sim3')
dataForbesAll = merge(dataForbesAll, dataForbes, by.x='trackName', by.y='trackName')

dataForbes = read.table(forbes_input4, head=F, sep='\t')[,c(2,3)]
colnames(dataForbes) = c('trackName', 'sim4')
dataForbesAll = merge(dataForbesAll, dataForbes, by.x='trackName', by.y='trackName')

dataForbes = read.table(forbes_input5, head=F, sep='\t')[,c(2,3)]
colnames(dataForbes) = c('trackName', 'sim5')
dataForbesAll = merge(dataForbesAll, dataForbes, by.x='trackName', by.y='trackName')

dataForbes = read.table(forbes_input6, head=F, sep='\t')[,c(2,3)]
colnames(dataForbes) = c('trackName', 'sim6')
dataForbesAll = merge(dataForbesAll, dataForbes, by.x='trackName', by.y='trackName')

dataForbes = read.table(forbes_input7, head=F, sep='\t')[,c(2,3)]
colnames(dataForbes) = c('trackName', 'sim7')
dataForbesAll = merge(dataForbesAll, dataForbes, by.x='trackName', by.y='trackName')

dataForbes = read.table(forbes_input8, head=F, sep='\t')[,c(2,3)]
colnames(dataForbes) = c('trackName', 'sim8')
dataForbesAll = merge(dataForbesAll, dataForbes, by.x='trackName', by.y='trackName')

dataForbes = read.table(forbes_input9, head=F, sep='\t')[,c(2,3)]
colnames(dataForbes) = c('trackName', 'sim9')
dataForbesAll = merge(dataForbesAll, dataForbes, by.x='trackName', by.y='trackName')

###################
#### Jaccard
dataJaccardAll = read.table(jaccard_input0, head=F, sep='\t')[,c(2,3)]
colnames(dataJaccardAll) = c('trackName', 'simAll')

dataJaccard = read.table(jaccard_input1, head=F, sep='\t')[,c(2,3)]
colnames(dataJaccard) = c('trackName', 'sim1')
dataJaccardAll = merge(dataJaccardAll, dataJaccard, by.x='trackName', by.y='trackName')

dataJaccard = read.table(jaccard_input2, head=F, sep='\t')[,c(2,3)]
colnames(dataJaccard) = c('trackName', 'sim2')
dataJaccardAll = merge(dataJaccardAll, dataJaccard, by.x='trackName', by.y='trackName')

dataJaccard = read.table(jaccard_input3, head=F, sep='\t')[,c(2,3)]
colnames(dataJaccard) = c('trackName', 'sim3')
dataJaccardAll = merge(dataJaccardAll, dataJaccard, by.x='trackName', by.y='trackName')

dataJaccard = read.table(jaccard_input4, head=F, sep='\t')[,c(2,3)]
colnames(dataJaccard) = c('trackName', 'sim4')
dataJaccardAll = merge(dataJaccardAll, dataJaccard, by.x='trackName', by.y='trackName')

dataJaccard = read.table(jaccard_input5, head=F, sep='\t')[,c(2,3)]
colnames(dataJaccard) = c('trackName', 'sim5')
dataJaccardAll = merge(dataJaccardAll, dataJaccard, by.x='trackName', by.y='trackName')

dataJaccard = read.table(jaccard_input6, head=F, sep='\t')[,c(2,3)]
colnames(dataJaccard) = c('trackName', 'sim6')
dataJaccardAll = merge(dataJaccardAll, dataJaccard, by.x='trackName', by.y='trackName')

dataJaccard = read.table(jaccard_input7, head=F, sep='\t')[,c(2,3)]
colnames(dataJaccard) = c('trackName', 'sim7')
dataJaccardAll = merge(dataJaccardAll, dataJaccard, by.x='trackName', by.y='trackName')

dataJaccard = read.table(jaccard_input8, head=F, sep='\t')[,c(2,3)]
colnames(dataJaccard) = c('trackName', 'sim8')
dataJaccardAll = merge(dataJaccardAll, dataJaccard, by.x='trackName', by.y='trackName')

dataJaccard = read.table(jaccard_input9, head=F, sep='\t')[,c(2,3)]
colnames(dataJaccard) = c('trackName', 'sim9')
dataJaccardAll = merge(dataJaccardAll, dataJaccard, by.x='trackName', by.y='trackName')

##### Plot
dataForbesAll[dataForbesAll==0] = NA
dataJaccardAll[dataJaccardAll==0] = NA


###Loading dataset for grouping the data by cell type (download the file from GitHub)
groups = read.csv('cellTypeGroups.csv', sep=';')
dataJaccardAll = merge(dataJaccardAll, groups, by.x='trackName', by.y='trackName')
dataJaccardAll$groupColor = 'grey'
dataJaccardAll$groupColor[dataJaccardAll$PrimaryImmuneCells=='x'] = 'red'
dataJaccardAll$groupColor[dataJaccardAll$Thymus=='x'] = 'orange'
dataJaccardAll$groupColor[dataJaccardAll$Brain=='x'] = 'yellow'

dataForbesAll = merge(dataForbesAll, groups, by.x='trackName', by.y='trackName')
dataForbesAll$groupColor = 'grey'
dataForbesAll$groupColor[dataForbesAll$PrimaryImmuneCells=='x'] = 'red'
dataForbesAll$groupColor[dataForbesAll$Thymus=='x'] = 'orange'
dataForbesAll$groupColor[dataForbesAll$Brain=='x'] = 'yellow'

dataForbesAllForPlot = t(dataForbesAll)
dataJaccardAllForPlot = t(dataJaccardAll)


#### Plot Figure 5
pdf('figure5.pdf', width = 10, height = 12)
par(mfrow=c(2,1))

matplot(dataForbesAllForPlot[2:11,], col=dataForbesAll$groupColor, axes=F, type='l',pch = 16, xlab="GWAS P-value threshold", ylab='Fold enrichment of SNPs in DHSs', ylim=c(0,16))
axis(1, 1:10, c('1.0', '0.1', '0.01', '0.001', '0.0001', '10^-5', '10^-7', '10^-10', '10^-15', '10^-20'))
axis(2, seq(0,16,2), seq(0,16,2))

matplot(dataJaccardAllForPlot[2:11,], col=dataJaccardAll$groupColor, axes=F, type='l', pch = 16, xlab="GWAS P-value threshold", ylab='Jaccard similarity of SNPs in DHSs', ylim=c(0,0.0005))
axis(1, 1:10, c('1.0', '0.1', '0.01', '0.001', '0.0001', '10^-5', '10^-7', '10^-10', '10^-15', '10^-20'))
axis(2, seq(0,0.0005, 0.00005), seq(0,0.0005, 0.00005))

legend(6, 0.0004, c('immune cells', 'thymus', 'brain', 'other'), col=c('red', 'orange', 'yellow', 'grey'), pch=1, box.col = 'white', lty=1)

dev.off()











