## package names
packages <- c("data.table", "readxl", "knitr", "IlluminaHumanMethylation450kanno.ilmn12.hg19")

## install cran packages
install.packages(head(packages, -1))

## install bioconductor package
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

## check that all were installed and can be loaded
sapply(packages, require, character.only=T)