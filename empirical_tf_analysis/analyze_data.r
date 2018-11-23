library(nlme)
.get.coefficient.name <- function(name) paste("group", name, sep="")

analyze.scores.adv <- function(data, name) {
    data$score.rank <- rank(data$score)
    data$size.sum.rank <- rank(data$size)
    data$size.prod.rank <- rank(data$logsize)
    model <- lm(score.rank ~ size.sum.rank + size.prod.rank + group, data = data)
    # model <- lm(score.rank ~ size.prod.rank + group, data = data)
    model.mixed <- lme(score.rank ~ size.sum.rank + size.prod.rank,random=~1|group, data = data)
    # model.mixed <- lme(score.rank ~ size.prod.rank,random=~1|group, data = data)
    print(summary(model.mixed))
    data$offsets <- model$coefficients[sapply(data$group, .get.coefficient.name)]
    data$offsets[is.na(data$offsets)] = 0
    group.sd <- diag(sqrt(getVarCov(model.mixed)))
    correlation <- group.sd/(group.sd+model.mixed$sigma)
    print("Correlation:")
    print(correlation)
    prod.res <- data$score.rank - data$offsets-data$size.sum.rank*model$coefficients["size.sum.rank"]
    plot(data$size.prod.rank, prod.res)
    abline(model$coefficients["(Intercept)"], model$coefficients["size.prod.rank"], col="red")
    sum.res <- data$score.rank-data$offsets-data$size.prod.rank*model$coefficients["size.prod.rank"]
    plot(data$size.sum.rank, sum.res)
    abline(model$coefficients["(Intercept)"], model$coefficients["size.sum.rank"], col="red")
    return(data.frame(x.prod=data$size.prod.rank, y.prod=prod.res, x.sum=data$size.sum.rank, y.sum=sum.res))
}

analyze.scores.predvres <- function(data, name) {
    data$score.rank <- rank(data$score)
    data$size.sum.rank <- rank(data$size)
    data$size.prod.rank <- rank(data$logsize)
    model <- lm(score.rank ~ size.sum.rank + size.prod.rank + group, data = data)
    model.mixed <- lme(score.rank ~ size.sum.rank + size.prod.rank,random=~1|group, data = data)
    print(summary(model.mixed))
    data$offsets <- model$coefficients[sapply(data$group, .get.coefficient.name)]
    data$offsets[is.na(data$offsets)] = 0
    plot(data$size.prod.rank*model$coefficients["size.prod.rank"]+
         data$size.sum.rank*model$coefficients["size.sum.rank"]), data$score.rank-data$offsets)
}

analyze.scores.raw <- function(data, name) {
    data$score.rank <- data$score # rank(data$score)
    data$size.sum.rank <- data$size
    data$size.prod.rank <- data$logsize
    model.prod <- lm(score.rank ~ size.prod.rank, data=data)
    model.sum <- lm(score.rank ~ size.sum.rank, data=data)
    plot(data$size.prod.rank, data$score.rank)
    abline(model.prod$coefficients["(Intercept)"], model.prod$coefficients["size.prod.rank"], col="red")
    plot(data$size.sum.rank, data$score.rank)
    abline(model.sum$coefficients["(Intercept)"], model.sum$coefficients["size.sum.rank"], col="red")
    return(data.frame(x.prod=data$size.prod.rank, x.sum=data$size.sum.rank, y=data$score.rank))
}

analyze.scores <- function(data, name) {
    data$score.rank <- rank(data$score) # rank(data$score)
    data$size.sum.rank <- rank(data$size)
    data$size.prod.rank <- rank(data$logsize)
    model.prod <- lm(score.rank ~ size.prod.rank, data=data)
    model.sum <- lm(score.rank ~ size.sum.rank, data=data)
    plot(data$size.prod.rank, data$score.rank)
    abline(model.prod$coefficients["(Intercept)"], model.prod$coefficients["size.prod.rank"], col="red")
    plot(data$size.sum.rank, data$score.rank)
    abline(model.sum$coefficients["(Intercept)"], model.sum$coefficients["size.sum.rank"], col="red")
    return(data.frame(x.prod=data$size.prod.rank, x.sum=data$size.sum.rank, y=data$score.rank))
}

analyze.scores.basic <- function(data) {
    data$score.rank <- rank(data$score)
    data$size.sum.rank <- rank(data$size)
    data$size.prod.rank <- rank(data$logsize)
    model <- lm(score.rank ~ size.prod.rank + group, data = data)
    # model <- lm(score.rank ~ size.prod.rank + group, data = data)
    model.mixed <- lme(score.rank ~ size.prod.rank,random=~1|group, data = data)
    # model.mixed <- lme(score.rank ~ size.prod.rank,random=~1|group, data = data)
    print(summary(model.mixed))
    data$offsets <- model$coefficients[sapply(data$group, .get.coefficient.name)]
    data$offsets[is.na(data$offsets)] = 0
    group.sd <- diag(sqrt(getVarCov(model.mixed)))
    correlation <- group.sd/(group.sd+model.mixed$sigma)
    print("Correlation:")
    print(correlation)
    plot(data$size.prod.rank,
         data$score.rank - data$offsets)
    abline(model$coefficients["(Intercept)"], model$coefficients["size.prod.rank"], col="red")
    return(data.frame(x=data$size.prod.rank, y=data$score.rank - data$offsets))
}


main.basic <- function(jaccard, forbes, tetra){
    png("basic.png")
    par(mfrow=c(3,1))
    write.csv(analyze.scores.basic(jaccard), "jaccard_basic_points.csv", row.names=F)
    write.csv(analyze.scores.basic(forbes) , "forbes_basic_points.csv", row.names=F)
    write.csv(analyze.scores.basic(tetra)  , "tetra_basic_points.csv", row.names=F)
    dev.off()
}

main.raw <- function(jaccard, forbes, tetra){
    png("raw.png")
    par(mfrow=c(3, 2))
    write.csv(analyze.scores.raw(jaccard), file="jaccard_raw_points.csv", row.names=F)
    write.csv(analyze.scores.raw(forbes), file="forbes_raw_points.csv", row.names=F)
    write.csv(analyze.scores.raw(tetra), file="tetra_raw_points.csv", row.names=F)
    dev.off()
}

main <- function(jaccard, forbes, tetra){
    png("rawrank.png")
    par(mfrow=c(3, 2))
    write.csv(analyze.scores(jaccard), file="jaccard_rawrank_points.csv", row.names=F)
    write.csv(analyze.scores(forbes), file="forbes_rawrank_points.csv", row.names=F)
    write.csv(analyze.scores(tetra), file="tetra_rawrank_points.csv", row.names=F)
    dev.off()
}

main.predvres <- function(jaccard, forbes, tetra){
    # png("predvres.png")
    par(mfrow=c(3, 1))
    analyze.scores.predvres(jaccard)
    analyze.scores.predvres(forbes)
    analyze.scores.predvres(tetra)
    # dev.off()
}


    # hlines <- rs-mean(rs)
    # hlines <- sort(hlines)
    # print(length(hlines))
    # print(seq(1, length(hlines), 10))
    # print(hlines[seq(1, length(hlines), 10)])
    # abline(model$coefficients["(Intercept)"], model$coefficients["size.prod.rank"], col="red")
    # for (line in hlines[seq(1, length(hlines), 10)]) abline(h=line, col="green")
    # abline(0, model$coefficients["size.prod.rank"], col="red")
    # abline(model$coefficients["(Intercept)"]+group.sd, model$coefficients["size.prod.rank"], col="green")
    # abline(model$coefficients["(Intercept)"]-group.sd, model$coefficients["size.prod.rank"], col="green")
    # plot(offsets-min(offsets), score.rank)
    # abline(model$coefficients["(Intercept)"]+group.sd, model$coefficients["size.sum.rank"], col="green")
    # abline(model$coefficients["(Intercept)"]-group.sd, model$coefficients["size.sum.rank"], col="green")
