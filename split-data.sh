#!/bin/bash 

source variables.txt

mkdir -p $PROTEIN_DIR
mkdir -p $METHYLATION_DIR

Rscript split-data.r \
    $FULL_DIR/clinical-clean.txt \
    $FULL_DIR/methylation-clean.txt \
    $FULL_DIR/protein-clean.txt \
    $PROTEIN_DIR \
    $METHYLATION_DIR


