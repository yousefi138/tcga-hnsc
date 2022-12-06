args <- c("hnsc-project/data/episcores-files/Predictors_Shiny_by_Groups.csv",
			"hnsc-project/data/methylation-dataset")

episcore.weights.filename <- args[1]
methylation.dir <- args[2]
methylation.file <- file.path(methylation.dir, "methylation.txt")

## Start to Process Files 

my.read.table <- function(filename, ...) {
	require(data.table)
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
data <- my.read.table(methylation.file)
	index <- grep("Hybrid", colnames(data))
	rownames(data) <- data[,index[1]]
	data <- as.matrix(data[, -index])


## check number of rows missing per sample
miss <- apply(data, 2, function(i) table(is.na(i)), simplify=F)
miss.df <- as.data.frame(do.call(rbind, miss))
summary(miss.df$"TRUE")

## not sure why these observations are missing ~90k observations
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  89519   89566   89645   89790   89817   95558


library(meffonym)
library(tidyverse)

models <- meffonym.models(full=T) %>%
            filter(str_detect(filename, "^episcores"))

proteins <- models %>% pull(name)


episcores <- sapply(
        proteins,
        function(model) {
            cat(date(), model, " ")
            ret <- meffonym.score(data, model)
            cat(" used ", length(ret$sites), "/", length(ret$vars), "sites\n")
            ret$score
        })

## issue arises when some cpgs are missing for all samples... can we just set these to zero?
# ue Dec  6 21:09:18 2022 CCL18   used  119 / 120 sites
# ue Dec  6 21:09:18 2022 CCL21   used  98 / 99 sites
# ue Dec  6 21:09:18 2022 CCL22  Error in impute.mean(x, 1, na.rm = T) : all(!is.na(x)) is not TRUE

