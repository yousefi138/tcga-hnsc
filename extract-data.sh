#!/bin/bash

source variables.txt

mkdir -p $FULL_DIR

#############################
## bash function for extracting a specific 
## file from a given tar.gz file. 
## Usage: extract_file TAR_FILE EXTRACT_FILE NEW_FILE
extract_file () {
  if [ $# != 3 ]; then
    echo "Usage: extract_file TAR_FILE EXTRACT_FILE NEW_FILE"
    exit 1
  fi
  if [ -f "$3" ]; then
    FILE=`basename $3`
    echo "file $FILE exists, skipping"
    return
  fi
  echo "extracting $2 => $3"
  TAR_FILE=$1
  EXTRACT_FILE=`tar -ztf $TAR_FILE | grep -e $2`
  NEW_FILE=$3
  tar -zxvf $TAR_FILE $EXTRACT_FILE
  mv $EXTRACT_FILE $NEW_FILE
  rm -rf `dirname $EXTRACT_FILE`
}

########################
## protein
extract_file \
 $FILE_DIR/*_protein_normalization__data.Level_3.*.tar.gz \
 data.txt \
 $FULL_DIR/protein.txt

if [ ! -f "$FULL_DIR/protein-clean.txt" ]; then
    echo "Cleaning protein data"
    sed 2d $FULL_DIR/protein.txt > $FULL_DIR/protein-clean.txt
fi

########################
## methylation
extract_file \
 $FILE_DIR/*_HNSC.*_humanmethylation450_*_data.Level_3.*.tar.gz \
 data.txt \
 $FULL_DIR/methylation.txt

if [ ! -f $FULL_DIR/methylation-clean.txt ]; then
    echo "Cleaning up methylation data"
    awk -F'\t' '{printf "%s\t",$1; for(i=2;i<=NF;i+=4){printf "%s\t",$i;} print ""}'  \
	$FULL_DIR/methylation.txt \
	| sed 2d \
	> $FULL_DIR/methylation-clean.txt
fi

#######################
## clinical
extract_file \
 $FILE_DIR/gdac.broadinstitute.org_HNSC.Merge_Clinical.Level_1.2016012800.0.0.tar.gz \
 HNSC.clin.merged.txt \
 $FULL_DIR/clinical.txt

if [ ! -f "$FULL_DIR/clinical-clean.txt" ]; then
    echo "Cleaning clinical data"
    Rscript clean-clinical.r \
	$FULL_DIR/clinical.txt \
	$PAN_CANCER_CLINICAL \
	$FULL_DIR/clinical-clean.txt
fi
