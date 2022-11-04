#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

protein.dir <- args[1]
methylation.dir <- args[2]

cat("annotate-data.r", paste(args,collapse=" "), "\n")

############################
## load all data types other than methylation
source("my-read-table.r")
protein <- my.read.table(file.path(protein.dir, "protein.txt"))


############################
## methylation data annotation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(list="IlluminaHumanMethylation450kanno.ilmn12.hg19")
data(Locations)
data(Other)
meth.annotation <- data.frame(
    as.data.frame(Locations, stringsAsFactors=F),
    as.data.frame(Other, stringsAsFactors=F),
    stringsAsFactors=F)
cols <- c(
    chr="chr",
    pos="pos",
    gene="UCSC_RefGene_Name",
    region="UCSC_RefGene_Group")
meth.annotation <- meth.annotation[,cols]
colnames(meth.annotation) <- names(cols)
meth.annotation$site <- rownames(meth.annotation)

source("expand-data-frame.r")

## The 'gene' and 'region' columns are semi-colon delimited lists
## of gene/region pairs. Expand meth.annotation
## so that there is one row for each element
## e.g. a row like the following:
##  "chrY 9363356 TSPY;FAM197Y2 Body;TSS1500 cg00050873"
## becomes two rows:
##  "chrY 9363356 TSPY     Body    cg00050873"
##  "chrY 9363356 FAM197Y2 TSS1500 cg00050873"
meth.annotation$gene <- strsplit(meth.annotation$gene, ";")
meth.annotation$region <- strsplit(meth.annotation$region, ";")
regions <- expand.data.frame(meth.annotation, "region")$region
meth.annotation <- expand.data.frame(meth.annotation, "gene")
meth.annotation$region <- regions

write.table(meth.annotation, file=file.path(methylation.dir, "methylation-annotation.txt"), row.names=F, sep="\t")


##############
## protein identifiers to gene symbols
protein.annotation <- read.table("protein-annotation.txt",header=T,sep="\t")

write.table(protein.annotation,file=file.path(protein.dir, "protein-annotation.txt"), row.names=F, sep="\t")

