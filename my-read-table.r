library(data.table)
my.read.table <- function(filename, ...) {
    cat("reading", basename(filename), "... ")
    ## read in tab-delimited spreadsheet
    x <- fread(
        filename,
        header=T,
        stringsAsFactors=F,
        sep="\t",
        ...)
    ## remove any duplicate rows (identified by the first column)
    x <- x[match(unique(x[[1]]), x[[1]]),]
    ## make the first column the rownames of the data frame
    x <- data.frame(x,row.names=1,stringsAsFactors=F)
    cat(nrow(x), "x", ncol(x), "\n")
    x
}
