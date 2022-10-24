#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

training.dir <- args[1]
testing.dir <- args[2]

cat("run-test.r", paste(args,collapse=" "), "\n")


source("my-read-table.r")

##############################################
## training

## list all files in the training dataset
filenames <- list.files(training.dir, pattern="txt$", full.names=T)
filenames <- filenames[!grepl("(annotation|targets)",filenames)]

## for this example we'll omit CNV
filenames <- filenames[!grepl("cnv", filenames)]

## load the data files into a list
dat <- lapply(filenames, my.read.table)

## name the items of the list by the filename
names(dat) <- sub(".txt", "", basename(filenames))

## remove the clinical data from the list,
## it's a bit different than the other data files
clinical <- dat$clinical
dat$clinical <- NULL

## outcome if progression-free interval
outcome.var <- "pfi"

## clinical variables that might be informative about outcome
clinical.vars <- c("age.at.diagnosis","estrogen.receptor.status",
                   "progesterone.receptor.status",
                   "lymphocyte.infiltration","necrosis.percent")

## remove features with no variance or all missing values
for (name in names(dat)) {
    feature.var <- apply(dat[[name]],1,var,na.rm=T)
    dat[[name]] <- dat[[name]][which(feature.var > 2e-16),]
}
    

## identify top univariate predictors in each data type
library(limma)
univariate.predictors <- sapply(names(dat), function(name) {
    ## for each data type ...
    cat(date(), "testing", name, "...\n")
    ## prepare to test, for each feature, feature ~ outcome
    outcome <- clinical$pfi[match(colnames(dat[[name]]),rownames(clinical))]
    design <- model.matrix(~outcome)
    ## fit linear model for each feature
    fit <- lmFit(dat[[name]], design)
    ## calculate p-values
    fit <- eBayes(fit)
    ## identify the top 100 associations
    idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
    ## return the names of the features with the top 25 associations
    rownames(dat[[name]])[idx]
})
names(univariate.predictors) <- names(dat)

## fit an elastic net on the samples with complete data
## including top univeriate predictors

## R package providing elastic net functionality
library(glmnet) 

## convert clinical data to a numeric matrix
##  (i.e. create dummy variables for categorical variables)
clinical <- model.matrix(~0+., clinical[,c(outcome.var, clinical.vars)])

## identify participants with data for all data types
common.ids <- rownames(clinical)
for (i in 1:length(dat))
    common.ids <- intersect(common.ids, colnames(dat[[i]]))

## construct a dataset including individuals with data for all data types
## and limit to features with strong univariate associations
univariate.dat <- sapply(names(dat), function(name) {
    t(dat[[name]][univariate.predictors[,name],common.ids])
}, simplify=F)
univariate.dat$clinical <- clinical[common.ids,setdiff(colnames(clinical),outcome.var)]

## merge data for each data type into a single matrix
univariate.dat <- do.call(cbind, univariate.dat)

## impute missing values with the median value for the feature
idx <- which(is.na(univariate.dat),arr.ind=T)
median.values <- apply(univariate.dat,2,median,na.rm=T)
univariate.dat[idx] <- median.values[idx[,2]]

## outcome variable
outcome <- clinical[common.ids,outcome.var]


## standardize features
feature.var <- apply(univariate.dat,2,var)
univariate.dat <- univariate.dat[,feature.var > 2e-16]
univariate.dat <- scale(univariate.dat)

## fit elastic net model
cvfit <- cv.glmnet(
    univariate.dat,
    outcome,
    type.measure = "class",
    nfolds = 20,
    family="binomial") ## the outcome variable is binary


## apply model in the training dataset
train.predicts <- predict(
    cvfit,
    newx=univariate.dat,
    s="lambda.min",
    type="response")

## check it's performance
## area under the curve
library(pROC)
auc(outcome, as.vector(train.predicts))
## Area under the curve: 0.9487

###############################################
## testing

## load the test dataset
filenames <- list.files(testing.dir, full.names=T)
filenames <- filenames[!grepl("cnv", filenames)] ## ignoring cnv
test.dat <- lapply(filenames, my.read.table)
names(test.dat) <- sub(".txt", "", basename(filenames))

## extract the clinical data
test.clinical <- test.dat$clinical
test.dat$clinical <- NULL

## convert clinical data to a numeric matrix
## (this mainly means replacing categorical variables with dummy variables)
test.clinical <- model.matrix(~0+., test.clinical[,c(outcome.var, clinical.vars)])

## identify participants with data for all data types
common.ids <- rownames(test.clinical)
for (i in 1:length(test.dat))
    common.ids <- intersect(common.ids, colnames(test.dat[[i]]))

## restrict test data to univariate predictor features
## identified in training and to individuals data for all data types
test.dat <- sapply(names(test.dat), function(name) {
    t(test.dat[[name]][univariate.predictors[,name],common.ids])
}, simplify=F)
test.dat$clinical <- test.clinical[common.ids,setdiff(colnames(test.clinical),outcome.var)]

## merge data types into a single data matrix
test.dat <- do.call(cbind, test.dat)

## standardize features
feature.medians <- apply(test.dat,2,median,na.rm=T)
features.var <- apply(test.dat, 2, var, na.rm=T)
test.dat <- scale(test.dat)

## impute missing values with the median value for the feature
idx <- which(is.na(test.dat),arr.ind=T)
test.dat[idx] <- feature.medians[idx[,2]]

## apply trained model in the test dataset
test.predicts <- predict(
    cvfit,
    newx = test.dat[,colnames(univariate.dat)],
    s = "lambda.min",
    type = "response")

## evaluate performance

## outcome to predict
test.outcome <- test.clinical[common.ids,outcome.var]

## area under the curve
auc(test.outcome,as.vector(test.predicts))
## Area under the curve: 0.5639
## ... hopefully you can do better!


## check agreement about univariate predictors between
## training and testing
tfit <- lmFit(t(test.dat), model.matrix(~test.outcome))
tfit <- eBayes(tfit)

fit <- lmFit(t(univariate.dat), model.matrix(~outcome))
fit <- eBayes(fit)

cor(tfit$coef[rownames(fit$coef),"test.outcome"], fit$coef[,"outcome"])
## [1] 0.2417875
cor(tfit$t[rownames(fit$t),"test.outcome"], fit$t[,"outcome"])
## [1] 0.2421426

