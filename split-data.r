#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

filenames <- head(args,-2)
names(filenames) <- c(
    "clinical",
    "methylation",
    "protein")

protein.dir <- args[length(args)-1]
methylation.dir <- args[length(args)]

dir.create(protein.dir, recursive=T)
dir.create(c, recursive=T)

cat("split-data.r", paste(args,collapse=" "), "\n")

source("extract-participant.r")

library(data.table)
my.read.table <- function(filename, ...) {
    cat("reading", basename(filename), "... ")
    x <- fread(
        filename,
        header=T,
        stringsAsFactors=F,
        sep="\t",
        ...)
    cat(nrow(x), "x", ncol(x), "\n")
    as.data.frame(x,stringsAsFactors=F)
}

## load datasets
datasets <- lapply(filenames, my.read.table)

bycol <- setdiff(names(datasets), c("clinical"))

## extract participant names from sample ids for each dataset
for (name in bycol) {
    samples <- colnames(datasets[[name]])
    colnames(datasets[[name]]) <- extract.participant(samples)
}
datasets$clinical$participant <- extract.participant(datasets$clinical$participant)

        
## restrict the datasets to a set of participants
restrict.datasets <- function(datasets, participants) {
    c(list(
        clinical=datasets$clinical[datasets$clinical$participant %in% participants,]),
      sapply(
          datasets[which(names(datasets) %in% bycol)],
          function(dat) {
              idx <- c(1,which(colnames(dat) %in% participants))
              dat[,idx]
          }, simplify=F))
}

## restrict to those who have clinical data
datasets <- restrict.datasets(datasets, datasets$clinical$participant)

## create protein and methylation dataset subsets

## protein
protein.dataset <- restrict.datasets(datasets, colnames(datasets$protein))[c("clinical", "protein")]

## methylation
methylation.dataset <- restrict.datasets(datasets, colnames(datasets$methylation))[c("clinical", "methylation")]



library(knitr)
describe.subset <- function(datasets) {
    participants <- lapply(datasets, colnames)
    participants$clinical <- datasets$clinical$participant
    int <- sapply(
        participants,
        function(a) {
            sapply(
                participants,
                function(b) length(intersect(a,b)))
        })            
    colMeans(int)
}
describe.subset(datasets)
describe.subset(protein.dataset)
describe.subset(methylation.dataset)
   ## clinical methylation     protein
   ## 422.6667    440.3333    212.3333

   ## clinical  protein
   ##  212      213

   ## clinical methylation
   ##  528         555    

complete.deaths <- function(datasets) {
    participants <- lapply(datasets, colnames)
    participants$clinical <- datasets$clinical$participant
    deaths <- datasets$clinical$participant[datasets$clinical$pfi==1]
    participants <- lapply(participants, function(x) intersect(x,deaths))
    table(table(unlist(participants)))
}

complete.deaths(datasets)
complete.deaths(protein.dataset)
complete.deaths(methylation.dataset)
 ##   2   3 
 ## 104  94 

 ##  2
 ## 94

 ##   2
 ## 198


my.write.table <- function(x, filename) {
    cat("saving", basename(filename), "...\n")
    write.table(x, file=filename, row.names=F, col.names=T, sep="\t")
}

save.subset <- function(datasets,dir) {
    for (name in names(datasets))
        my.write.table(datasets[[name]], file.path(dir,paste0(name,".txt")))
}

save.subset(protein.dataset, protein.dir)

save.subset(methylation.dataset, methylation.dir)
