#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# VCF_ALL_Beagle_pipeline.sh
# --------------------------
# SNP pipeline for Beagle (ALL population).
# (LOCAL SCRIPT)
#***************************************************************************

# The user must specify the path to the repository and GATK
REPOSITORY_PATH="/path/to/repository"
GATK_PATH="/path/to/gatk"

python ./local/2_FilterGenotypes.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population ALL -version Tcas3.30 -suffix Beagle
python ./local/3_VariantsToTable.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population ALL -version Tcas3.30 -suffix Beagle -annotation

