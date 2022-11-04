#!/bin/bash

source variables.txt

Rscript annotate-data.r $PROTEIN_DIR $METHYLATION_DIR

