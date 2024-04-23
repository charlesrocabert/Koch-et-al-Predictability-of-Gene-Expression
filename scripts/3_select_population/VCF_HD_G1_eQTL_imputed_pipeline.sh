#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# VCF_HD_G1_eQTL_imputed_pipeline.sh
# ----------------------------------
# SNP pipeline for eQTL (HD G1 imputed population).
# (LOCAL SCRIPT)
#***************************************************************************

# The user must specify the path to the repository and GATK
REPOSITORY_PATH="/path/to/repository"
GATK_PATH="/path/to/gatk"

python ./local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix imputed -out-population HD_G1 -out-suffix imputed -version Tcas3.30 -task VCF
python ./local/2_FilterGenotypes.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population HD_G1 -version Tcas3.30 -imputed -suffix eQTL_imputed
python ./local/3_VariantsToTable.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population HD_G1 -version Tcas3.30 -suffix eQTL_imputed -annotation

