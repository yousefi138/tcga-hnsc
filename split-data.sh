#!/bin/bash 

source variables.txt

mkdir -p $TRAINING_DIR
mkdir -p $TESTING_DIR

Rscript split-data.r \
    $FULL_DIR/clinical-clean.txt \
    $FULL_DIR/mrna-clean.txt \
    $FULL_DIR/mirna.txt \
    $FULL_DIR/methylation-clean.txt \
    $FULL_DIR/protein-clean.txt \
    $FULL_DIR/mutations.txt \
    $FULL_DIR/cnv.txt \
    $TRAINING_DIR \
    $TESTING_DIR


