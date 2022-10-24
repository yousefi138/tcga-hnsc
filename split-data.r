#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

filenames <- head(args,-2)
names(filenames) <- c(
    "clinical",
    "mrna",
    "mirna",
    "methylation",
    "protein",
    "mutation",
    "cnv")

training.dir <- args[length(args)-1]
testing.dir <- args[length(args)]

dir.create(training.dir, recursive=T)
dir.create(testing.dir, recursive=T)

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

bycol <- setdiff(names(datasets), c("clinical","cnv"))

## extract participant names from sample ids for each dataset
for (name in bycol) {
    samples <- colnames(datasets[[name]])
    colnames(datasets[[name]]) <- extract.participant(samples)
}
datasets$clinical$participant <- extract.participant(datasets$clinical$participant)
datasets$cnv$Sample <- extract.participant(datasets$cnv$Sample)

        
## restrict the datasets to a set of participants
restrict.datasets <- function(datasets, participants) {
    c(list(
        clinical=datasets$clinical[datasets$clinical$participant %in% participants,],
        cnv=datasets$cnv[datasets$cnv$Sample %in% participants,]),
      sapply(
          datasets[which(names(datasets) %in% bycol)],
          function(dat) {
              idx <- c(1,which(colnames(dat) %in% participants))
              dat[,idx]
          }, simplify=F))
}

## restrict to those who have clinical data
datasets <- restrict.datasets(datasets, datasets$clinical$participant)


## create training and testing subsets

## select training
set.seed(20211112)
training.size <- floor(nrow(datasets$clinical)*2/3)
training.participants <- sample(datasets$clinical$participant,training.size)

## select testing
testing.participants <- setdiff(datasets$clinical$participant,training.participants)

## datasets
training.datasets <- restrict.datasets(datasets, training.participants)
testing.datasets <- restrict.datasets(datasets, testing.participants)


library(knitr)
describe.subset <- function(datasets) {
    participants <- lapply(datasets, colnames)
    participants$clinical <- datasets$clinical$participant
    participants$cnv <- unique(datasets$cnv$Sample)
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
describe.subset(training.datasets)
describe.subset(testing.datasets)
   ## clinical         cnv        mrna       mirna methylation     protein 
   ## 991.5714    991.5714   1023.1429    977.7143    772.0000    844.0000 
   ## mutation 
   ## 895.1429 
   ## clinical         cnv        mrna       mirna methylation     protein 
   ## 660.4286    660.4286    659.1429    652.4286    487.8571    550.1429 
   ## mutation 
   ## 599.5714 
   ## clinical         cnv        mrna       mirna methylation     protein 
   ## 331.1429    331.1429    330.0000    325.4286    254.5714    276.2857 
   ## mutation 
   ## 295.8571 


complete.deaths <- function(datasets) {
    participants <- lapply(datasets, colnames)
    participants$clinical <- datasets$clinical$participant
    participants$cnv <- unique(datasets$cnv$Sample)
    deaths <- datasets$clinical$participant[datasets$clinical$pfi==1]
    participants <- lapply(participants, function(x) intersect(x,deaths))
    table(table(unlist(participants)))
}

complete.deaths(datasets)
complete.deaths(training.datasets)
complete.deaths(testing.datasets)
 ## 5  6  7 
 ## 9 62 73 

 ## 5  6  7 
 ## 4 40 46 

 ## 5  6  7 
 ## 5 22 27 


my.write.table <- function(x, filename) {
    cat("saving", basename(filename), "...\n")
    write.table(x, file=filename, row.names=F, col.names=T, sep="\t")
}

save.subset <- function(datasets,dir) {
    for (name in names(datasets))
        my.write.table(datasets[[name]], file.path(dir,paste0(name,".txt")))
}

save.subset(training.datasets, training.dir)

save.subset(testing.datasets, testing.dir)
