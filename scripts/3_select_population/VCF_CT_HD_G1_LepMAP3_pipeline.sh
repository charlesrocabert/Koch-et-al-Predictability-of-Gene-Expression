#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# VCF_CT_HD_G1_LepMAP3_pipeline.sh
# --------------------------------
# SNP pipeline for Lep-MAP3 (CT/HD G1 population).
# (LOCAL SCRIPT)
#***************************************************************************

# The user must specify the path to the repository and GATK
REPOSITORY_PATH="/path/to/repository"
GATK_PATH="/path/to/gatk"

python ./local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix raw_SNP -out-population CT_HD_G1 -out-suffix raw_SNP -version Tcas3.30 -task VCF
python ./local/2_FilterGenotypes.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population CT_HD_G1 -version Tcas3.30 -suffix LepMAP3
python ./local/3_VariantsToTable.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population CT_HD_G1 -version Tcas3.30 -suffix LepMAP3 -annotation

