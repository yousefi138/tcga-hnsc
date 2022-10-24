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
### mrna
#extract_file \
# $FILE_DIR/*mRNAseq_Preprocess.Level_3.*.tar.gz \
# RSEM_Z_Score.txt \
# $FULL_DIR/mrna.txt
#
#if [ ! -f  $FULL_DIR/mrna-clean.txt ]; then
#    echo "Cleaning mrna.txt ..."
#    perl -pe 's/\|[0-9]+//'  $FULL_DIR/mrna.txt \
#	| grep -ve "?" >  $FULL_DIR/mrna-clean.txt
#    
#fi
#
########################
### mirna 
#extract_file \
# $FILE_DIR/*HNSC.*illuminaga_*_miR_gene_expression__data.Level_3.*.tar.gz \
# miR_gene_expression__data.data.txt \
# $FULL_DIR/mirna1.txt
#
#extract_file \
# $FILE_DIR/*HNSC.*illuminahiseq_*_miR_gene_expression__data.Level_3.*.tar.gz \
# miR_gene_expression__data.data.txt \
# $FULL_DIR/mirna2.txt
#
#if [ ! -f "$FULL_DIR/mirna.txt" ]; then
#    echo "Cleaning mirna data"
#    Rscript clean-mirna.r \
#	$FULL_DIR/mirna1.txt \
#	$FULL_DIR/mirna2.txt \
#	$FULL_DIR/mirna.txt
#fi
#
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

########################
### cnv
#extract_file \
# $FILE_DIR/*_HNSC.Merge_*_cnv_hg19__seg.Level_3.*.tar.gz \
# seg.txt \
# $FULL_DIR/cnv.txt
#
########################
### mutations
#if [ -d "$FULL_DIR/mutations" ]; then
#    echo "directory mutations/ exists, skipping"
#else
#    gunzip -c $FILE_DIR/*_HNSC.Mutation_Packager_Oncotated_Calls.Level_3.*.tar.gz | tar xv
#    mv gdac.broad*  $FULL_DIR/mutations
#fi
#
#if [ -f "$FULL_DIR/mutations.txt" ]; then
#    echo "file mutations.txt exists, skipping"
#else
#    Rscript extract-mutations.r $FULL_DIR/mutations $FULL_DIR/mutations.txt
#fi
#    

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
