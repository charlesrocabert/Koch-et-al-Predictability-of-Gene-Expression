#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 6_merge_eQTLs_pipeline.sh
# -------------------------
# Merge eQTL associations.
# (LOCAL SCRIPT)
#***************************************************************************


# The user must specify the path to the repository and GATK
REPOSITORY_PATH="/path/to/repository"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1) Compute the correlation between phenotype and fitness #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Compute the correlation between phenotype and fitness"
Rscript ./scripts/8_eQTLs/local/ComputeCorrelationToFitness.R $REPOSITORY_PATH CT_G1 Tcas3.30 EXPRESSION
Rscript ./scripts/8_eQTLs/local/ComputeCorrelationToFitness.R $REPOSITORY_PATH HD_G1 Tcas3.30 EXPRESSION

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2) Extract the list of gene positions                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Extract the list of gene positions"
python ./scripts/8_eQTLs/local/ExtractGenePos.py -repository-path $REPOSITORY_PATH -version Tcas3.30

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3) Merge eQTLs data with annotations                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Merge eQTLs data with annotations and fitness (imputed data)"
Rscript ./scripts/8_eQTLs/local/MergeEQTLsDatasets.R $REPOSITORY_PATH CT_G1 Tcas3.30 _imputed EXPRESSION CT
Rscript ./scripts/8_eQTLs/local/MergeEQTLsDatasets.R $REPOSITORY_PATH HD_G1 Tcas3.30 _imputed EXPRESSION HD

