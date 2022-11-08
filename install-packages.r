## package names
cran <- c("data.table", "readxl", "knitr")
bioc <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"

installed <- sapply(c(cran,bioc), require, character.only=T)
installed

if (sum(installed) < length(installed)){

    index <- intersect(
                which(names(installed) %in% cran), 
                which(installed == FALSE))

        if(length(index) > 0){
            cat("Installing", names(installed[index]), "from CRAN", "\n")            
            install.packages(names(installed[index]))
        }

    index <- intersect(
                which(names(installed) %in% bioc), 
                which(installed == FALSE))

        if(length(index) > 0){
            if (!require("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
            cat("Installing", names(installed[index]), "from Bioconductor", "\n")            
            BiocManager::install(names(installed[index]))
        }

}

installed <- sapply(c(cran,bioc), require, character.only=T)
cat("Of the below packages required: \n",
    names(installed),
    "\n \n the following are available: \n", 
    names(installed[which(installed == T)]), "\n")

