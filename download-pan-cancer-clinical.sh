#!/bin/bash

source variables.txt

mkdir -p $PAN_CANCER_DIR

Rscript download-pan-cancer-clinical.r $PAN_CANCER_URL $PAN_CANCER_CLINICAL
