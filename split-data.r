#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

filenames <- head(args,-2)
names(filenames) <- c(
    "clinical",
    "methylation",
    "protein")

protein.dir <- args[length(args)-1]
methylation.dir <- args[length(args)]

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

## extract participant tissue information 
tissues <- data.frame(
                participant = extract.participant(colnames(datasets$methylation)),
                tissue = extract.tissue(colnames(datasets$methylation)),
                participant.tissue = paste(extract.participant(colnames(datasets$methylation)), 
                                        extract.tissue(colnames(datasets$methylation)), sep = "-")
    )
tissues <- subset(tissues, tissue!= "06")

## extract participant and tissue names from sample ids for each dataset

#protein
samples <- colnames(datasets[["protein"]])
colnames(datasets[["protein"]]) <- extract.participant(samples)

#methylation
samples <- colnames(datasets[["methylation"]])
colnames(datasets[["methylation"]]) <- paste(extract.participant(samples), 
                                        extract.tissue(samples), sep = "-")

datasets$clinical$participant <- extract.participant(datasets$clinical$participant)


# protein 
idx <- c(1,which(colnames(datasets$protein) %in% datasets$clinical$participant))
protein.dataset <- list(
                    clinical=datasets$clinical[datasets$clinical$participant %in% colnames(datasets$protein),],
                    protein = datasets$protein[,idx])


# methylation
datasets$clinical <- merge(datasets$clinical, tissues, by.x = "participant")
datasets$clinical$tumor.or.normal <- ifelse(as.numeric(datasets$clinical$tissue) < 9, "tumor", "normal")

idx <- c(1,which(colnames(datasets$methylation) %in% datasets$clinical$participant.tissue))
methylation.dataset <- list(
                    clinical=datasets$clinical[datasets$clinical$participant.tissue %in% colnames(datasets$methylation),],
                    methylation = datasets$methylation[,idx])



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
