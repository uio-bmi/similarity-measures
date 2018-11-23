library(nlme)

GetBeta <- function(title, data) {
    df <- data[data$title == title,]
    model <- lm(rank(score) ~ log(size2) + name2, data=df)
    return(model$coefficients["log(size2)"])
}

GetNaiveVar <- function(title, data) {
    df <- data[data$title == title,]
    model <- lm(rank(score) ~ name2, data=df)
    return(Wie
    N = length(model$coefficients)
    vs <- model$coefficients[2: N]
    m <- sum(vs)/N
    return(sqrt(sum((vs-m)^2)/(N-1)))
}

WierdVar <- function(coefs) {
    N = length(coefs)
    vs <- coefs[2: N]
    m <- sum(vs)/N
    return(sqrt(sum((vs-m)^2)/(N-1)))
}

GetNaiveVar <- function(title, data) {
    df <- data[data$title == title,]
    model <- lm(rank(score) ~ name2, data=df)
    return(WierdVar(model$coefficients)
}

GetAdjustedVar <- function(title, data) {
    df <- data[data$title == title,]
    model <- lm(rank(score) ~ rank(size2) + name2, data=df)
    return(WierdVar(model$coefficients[2:length(model$coefficients)]))
    #model.mixed <- lme(rank(score) ~ rank(size2),random=~1|name2, data = df)
    #return(diag(sqrt(getVarCov(model.mixed))))
}


GetSizes <- function(title, data) {
    df <- data[data$title == title,]
    return(mean(df$size))
}

Analysis <- function(data) {
    data = data[data$title != "MACS_ENCSR000DZM_GM12878_STAT1", ]
    titles <- unique(data$title)
    betas <-  sapply(titles, GetBeta, data)
    sizes <-  sapply(titles, GetSizes, data)
    plot(log(sizes), betas)
    naives <- sapply(titles, GetNaiveVar, data)
    adj <- sapply(titles, GetAdjustedVar, data)
    plot(naives, adj)
}
