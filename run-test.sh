#!/bin/bash

source variables.txt

/usr/bin/Rscript run-test.r $TRAINING_DIR $TESTING_DIR
