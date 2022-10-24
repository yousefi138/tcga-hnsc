#!/bin/bash

source variables.txt

FILE=$1

Rscript download-data.r $FILE $FILE_DIR
