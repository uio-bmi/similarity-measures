
rankDataFrame <- function(data){
    data$score.rank <- rank(data$score)
    data$size.sum.rank <- rank(data$size)
    data$size.prod.rank <- rank(data$logsize)
    return (data)
}
.get.coefficient.name <- function(name) paste("group", name, sep="")

getOffsets <- function(data, y, effects) {
    formula <- as.formula(
        paste(y, paste(effects, collapse="+"), sep ="~"))
    model <- lm(formula, data = data)
    for (effect in effects) {
        if (effect == "group"){
            data[paste(effect, "offset", sep=".")] <- coefficients(model)[sapply(data[effect], .get.coefficient.name)]
        }
        else {
            data[paste(effect, "offset", sep=".")] <- data[effect]*coefficients(model)[effect]
        }
    }
    return(data)
}

pad <- function(x, width) {
    r = (width-1)/2
    # start = x[1:r]-mean(x[(r-10):r])+mean(x[1:10])
    start = rep(mean(x[1:2]), r)
    l = length(x)
    end = rep(mean(x[(l-1):l]), r)
    # end = x[(l-r+1):l]-mean(x[(l-r):(l-r+10)])+mean(x[(l-10):l])
    return(c(start, x, end))
}

rollMean <- function(x, width) {
    c <- cumsum(x)
    l <- length(x)
    return((c[(width):l]-c[1:(l-width+1)])/width)
}

rollStd <- function(x, width) {
    means <- rollMean(x, width)
    c <- cumsum(x**2)
    l <- length(x)
    sums = (c[(width):l]-c[1:(l-width+1)])/(width-1)
    return(sqrt(sums-means**2))
}


analyze.scores.adv <- function(data) {
    data <- rankDataFrame(data)
    effects = c("group", "size.prod.rank", "size.sum.rank")
    data <- getOffsets(data, effects)
    plot(data$size.prod.rank, 
         data$score.rank-data$size.sum.rank.offset-data$group.offset)
    plot(data$size.sum.rank, data$score.rank-data$size.prod.rank.offset-data$group.offset)
    return(data)
}

simple.figure <- function(data) {
    data$simple.x <- rank(data$logsize)
    data$simple.y <- 11-log.rank(data$score)
    # plot(data$simple.x, data$simple.y)
    scatter.smooth(data$simple.x, data$simple.y, span=0.1)
    return(data)
}

complicated.figure <- function(data) {
    data <- getOffsets(data, "score", c("group"))
    data$group.offset[is.na(data$group.offset)] = 0
    data$complicated.x <- rank(data$logsize)
    data$complicated.y <- 11-log.rank(data$score-data$group.offset)

    # plot(data$complicated.x, data$complicated.y)
    scatter.smooth(data$complicated.x, data$complicated.y, span=0.1)
    return(data)
}
score.rank.mean.adv <- function(data, width) {
    data <- rankDataFrame(data)
    effects = c("group", "size.prod.rank")
    data <- getOffsets(data, "score.rank", effects)
    data$group.offset[is.na(data$group.offset)] = 0
    data <- data[order(data$size.prod.rank),]
    score.ranks = data$score.rank-data$group.offset
    local.mean.rank <- rollMean(score.ranks, width)
    sds = rollStd(score.ranks, width)
    y.max <- max(local.mean.rank+sds)
    y.min <- min(local.mean.rank-sds)
    x = data$size.prod.rank[(width):length(data$size.prod.rank)]
    plot(x, local.mean.rank, lty="dashed", ylim=c(0, 40000))
    lines(x, local.mean.rank+sds)
    lines(x, local.mean.rank-sds)
    return(list(x=x, y=local.mean.rank, upper=local.mean.rank+sds, lower=local.mean.rank-sds, point.x=data$size.prod.rank, point.y=data$score.rank-data$group.offset))
}


score.rank.mean <- function(data, width) {
    data <- rankDataFrame(data)
    data <- data[order(data$size.prod.rank),]
    local.mean.rank <- rollMean(data$score.rank, width)
    sds = rollStd(data$score.rank, width)
    y.max <- max(local.mean.rank+sds)
    y.min <- min(local.mean.rank-sds)
    print(length(local.mean.rank))
    print(length(data$size.prod.rank))
    x = data$size.prod.rank[(width+1):length(data$size.prod.rank)]
    plot(x, local.mean.rank, lty="dashed", ylim=c(y.min, y.max))
    lines(x, local.mean.rank+sds)
    lines(x, local.mean.rank-sds)
    return(data.frame(x=x, y=local.mean.rank, upper=local.mean.rank+sds, lower=local.mean.rank-sds))
}

runData <- function() {
    jaccard <- read.table("data/jaccard.csv", header=T)
    forbes <- read.table("data/forbes.csv", header=T)
    tetra <- read.table("data/tetra.csv", header=T)
    # png("logrank.png")
    par(mfrow=c(3, 2))

    data.jaccard <- simple.figure(jaccard)
    data.jaccard <- complicated.figure(data.jaccard)
    write.csv(data.jaccard, file="jaccard_logrank.csv", row.names=F)

    data.forbes <- simple.figure(forbes)
    data.forbes <- complicated.figure(data.forbes)
    write.csv(data.forbes, file="forbes_logrank.csv", row.names=F)


    data.tetra <- simple.figure(tetra)
    data.tetra <- complicated.figure(data.tetra)

    write.csv(data.tetra, file="tetra_logrank.csv", row.names=F)
    # dev.off()
}

summary.output = function(data, means, outfile) {
    cat(means$x, file=outfile, sep="\t")
    cat("\n", file=outfile, append=T)
    cat(means$y, file=outfile, append=T, sep="\t")
    cat("\n", file=outfile, append=T)
    cat(means$upper, file=outfile, append=T, sep="\t")
    cat("\n", file=outfile, append=T)
    cat(means$lower, file=outfile, append=T, sep="\t")
    cat("\n", file=outfile, append=T)
    cat(means$point.x, file=outfile, append=T, sep="\t")
    cat("\n", file=outfile, append=T)
    cat(means$point.y, file=outfile, append=T, sep="\t")
# 
#     cat(rank(data$logsize), file=outfile, append=T, sep="\t")
#     cat("\n", file=outfile, append=T)
#     cat(rank(data$score), file=outfile, append=T, sep="\t")
# 
}


runAvgs <- function() {
    jaccard <- read.table("data/jaccard.csv", header=T)
    forbes <- read.table("data/forbes.csv", header=T)
    tetra <- read.table("data/tetra.csv", header=T)
    png("localavg1000.png")
    par(mfrow=c(3, 1))
    data.jaccard = score.rank.mean.adv(jaccard, 1001)
    print(summary(data.jaccard))
    summary.output(jaccard, data.jaccard, "jaccard_localavg1000.tsv")
    # write.csv(data.jaccard, file="jaccard_localavg1000.csv", row.names=F)
    
    data.forbes = score.rank.mean.adv(forbes, 1001)
    summary.output(forbes, data.forbes, "forbes_localavg1000.tsv")
    # write.csv(data.forbes, file="forbes_localavg1000.csv", row.names=F)

    data.tetra = score.rank.mean.adv(tetra, 1001)
    summary.output(tetra, data.tetra, "tetra_localavg1000.tsv")
    #write.csv(data.tetra, file="tetra_localavg1000.csv", row.names=F)
    dev.off()
}




log.rank <- function(values) log(length(values)+1-rank(values))



