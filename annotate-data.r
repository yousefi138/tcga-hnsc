#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

data.dir <- args[1]

cat("annotate-data.r", paste(args,collapse=" "), "\n")

############################
## load all data types other than methylation
source("my-read-table.r")
mirna <- my.read.table(file.path(data.dir, "mirna.txt"))
mrna <- my.read.table(file.path(data.dir, "mrna.txt"))
protein <- my.read.table(file.path(data.dir, "protein.txt"))
mutation <- my.read.table(file.path(data.dir, "mutation.txt"))


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

write.table(meth.annotation, file=file.path(data.dir, "methylation-annotation.txt"), row.names=F, sep="\t")


#####################
## microRNA targets
url <- "http://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Targets_Context_Scores.default_predictions.txt.zip"
filename <- file.path(data.dir, "mirna-targets.txt.zip")
download.file(url, destfile=filename)
library(readr)
targets <- readr::read_tsv(filename)
targets <- as.data.frame(targets)
targets <- targets[,c("Gene Symbol","miRNA","context++ score")]
colnames(targets) <- c("gene","mirna","score")
targets <- targets[grepl("^hsa",targets$mirna),]
targets$mirna <- tolower(targets$mirna)

library("miRBaseConverter")
mirna.table <- miRBaseConverter::getMiRNATable(species="hsa")
idx <- match(targets$mirna, tolower(mirna.table$Mature1))
targets$precursor <- tolower(mirna.table$Precursor[idx])
is.missing <- is.na(targets$precursor)
idx <- match(targets$mirna[is.missing], tolower(mirna.table$Mature2))
targets$precursor[is.missing] <- tolower(mirna.table$Precursor[idx])

write.table(targets, file=file.path(data.dir, "mirna-targets.txt"), row.names=F, sep="\t")

##############
## protein identifiers to gene symbols
protein.annotation <- read.table("protein-annotation.txt",header=T,sep="\t")

write.table(protein.annotation,file=file.path(data.dir, "protein-annotation.txt"), row.names=F, sep="\t")

###############
## gene annotation

library(biomaRt)
mart <- useEnsembl(
        biomart="ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl",
        GRCh=37)
genes <- getBM(attributes = c("ensembl_gene_id",
                   "hgnc_symbol",
                   "chromosome_name",
                   "transcription_start_site",
                   "transcript_start",
                   "transcript_end",
                   "strand",
                   "gene_biotype"),
               mart = mart) ## 1 minute
genes <- genes[genes$gene_biotype=="protein_coding",]
genes <- genes[genes$chromosome_name %in% c(1:22,"X","Y"),]

genes <- genes[genes$hgnc_symbol %in% rownames(mrna),]
genes$ensembl_gene_id <- NULL
genes$gene_biotype <- NULL

colnames(genes) <- c("gene","chr","transcription_start_site","transcript_start","transcript_end","strand")

write.table(genes, file=file.path(data.dir, "gene-annotation.txt"),
            sep="\t", row.names=F)

##################
## check gene identifiers match 'enough' between data types

## methylation genes <-> mrna data
mean(meth.annotation$gene %in% rownames(mrna))
## 0.8360575

## mirna targets <-> mirna data
sum(rownames(mirna) %in% targets$precursor)/length(unique(targets$precursor))
## [1] 0.931

## mirna targets <-> mrna data
mean(targets$gene %in% rownames(mrna))
## [1] 0.943

## mrna data <-> protein data
mean(protein.annotation$gene %in% rownames(mrna))
## [1] 0.9912

## mutation data <-> mrna data
mean(rownames(mutation) %in% rownames(mrna))
## [1] 0.909

## gene annotation <-> mrna data
mean(rownames(mrna) %in% genes$gene)
## [1] 0.857
