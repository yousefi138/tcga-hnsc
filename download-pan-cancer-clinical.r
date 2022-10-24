#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

url <- args[1]
output.filename <- args[2]

dir.create(dirname(output.filename), recursive=T)

cat("download-pan-cancer-clinical.r", paste(args,collapse=" "), "\n")

xlsx.filename <- sub("txt$","xlsx",output.filename)

download.file(
    url,
    destfile=xlsx.filename)

library(readxl)

dat <- read_xlsx(xlsx.filename, sheet=1)

write.table(dat,file=output.filename,sep="\t",row.names=F,col.names=T)

