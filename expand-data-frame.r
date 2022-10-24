## expand a data frame so that a given column
## no longer has lists as elements
expand.data.frame <- function(v, column) {
    v[[column]][which(sapply(v[[column]], length) == 0)] <- NA
    x <- expand(v[[column]])
    w <- v[x$idx,]
    w[,column] <- x$value
    w
}

expand <- function(y) {
    idx <- rep(1:length(y), sapply(y, length))
    data.frame(idx=idx, value=unlist(y))
}

