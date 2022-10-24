#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

input.dir <- args[1]
output.filename <- args[2]

cat("extract-mutations.r", paste(args,collapse=" "), "\n")

source("extract-participant.r")

filenames <- list.files(input.dir, pattern="TCGA-.*.maf.txt", full.names=T)

## read a list of mutations for each sample
dat <- lapply(filenames, function(filename) {
    cat("reading", filename, "\n")
    read.table(
        filename,
        header=T,
        skip=3,
        sep="\t",
        quote="",
        comment.char="",
        stringsAsFactors=F)
})

dat <- do.call(rbind,dat)


dim(dat)
## [1] 90245   324

## remove silent mutations and common SNPs
dat.pruned <- dat[which(dat$Variant_Classification != "Silent"
                        & (dat$dbSNP_RS == "" | is.na(dat$dbSNP_RS))),]

dim(dat.pruned)
## [1] 62393   324

## list of genes with any mutations
symbols <- unique(dat.pruned$Hugo_Symbol)

## obtain participant ids
dat.pruned$participant <- extract.participant(dat.pruned$Tumor_Sample_Barcode)

## construct mutation profiles
## (rows=genes, columns=samples, values=number mutations in gene)
mutations <- sapply(unique(dat.pruned$participant), function(id) {
    mutated.genes <- dat.pruned$Hugo_Symbol[dat.pruned$participant == id]
    mutated.genes <- table(mutated.genes)
    profile <- rep(0,length(symbols))
    profile[match(names(mutated.genes), symbols)] <- mutated.genes
    profile
})
rownames(mutations) <- symbols

## how individuals with mutations per gene
gene.counts <- apply(mutations>0,1,sum)
quantile(gene.counts)
#  0%  25%  50%  75% 100% 
#   1    1    3    4  162 

mutations <- data.frame(gene=rownames(mutations), mutations, stringsAsFactors=F)
write.table(mutations, file=output.filename, row.names=F, col.names=T,sep="\t")

