#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# VCF_ALL_imputation_pipeline.sh
# ------------------------------
# Genotypes imputation pipeline (ALL population).
# (LOCAL SCRIPT)
#***************************************************************************

# The user must specify the path to the repository and Beagle
REPOSITORY_PATH="/path/to/repository"
BEAGLE_PATH="/path/to/beagle"

python ./local/3_ImputeGenotypes.py -repository-path $REPOSITORY_PATH -beagle $BEAGLE_PATH -population ALL -version Tcas3.30 -suffix Beagle

