#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

mirna.filenames <- c(args[1], args[2])
output.filename <- args[3]

cat("clean-mirna.r", paste(args,collapse=" "), "\n")

source("extract-participant.r")

dat <- lapply(mirna.filenames, function(filename)  {
    ## read first line (contains sample identifiers)
    samples <- readLines(filename,n=1)
    ## split into sample names (remove the first item, not a sample)
    samples <- strsplit(samples, "\t")[[1]][-1]
    ## extract participant ids
    samples <- extract.participant(samples)
    ## read the file (ignoring the first line, it just contains sample ids)
    dat <- read.table(filename,sep="\t",header=T,row.names=1,skip=1)
    ## there are multiple columns per sample, we just want reads/million
    is.count <- grepl("reads_per_million_miRNA_mapped",colnames(dat))
    ## extract just these
    dat <- dat[,is.count]
    ## replace column names with person identifiers
    colnames(dat) <- samples[is.count]    
    dat
})

## ensure rownames (gene names) match exactly between data files
stopifnot(
    all(sapply(
        dat,
        function(dati) identical(rownames(dat[[1]]), rownames(dati)))))

dat <- do.call(cbind, dat)

## 108 individuals have multiple miRNA profiles
table(table(colnames(dat)))
##   1   2   3 
## 971 105   3

## we retain the replicate with the highest correlation
## with the other profiles
r <- cor(dat)
r <- apply(r,2,median)
ids <- colnames(dat)
idx <- order(r,decreasing=T)
dat <- dat[,idx]
ids <- ids[idx]
dat <- dat[,match(unique(ids),ids)]

dat <- data.frame(mirna=rownames(dat), dat, stringsAsFactors=F)

write.table(dat, file=output.filename, row.names=F, col.names=T, sep="\t")

