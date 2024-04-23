#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# variant_detection_pipeline.sh
# -----------------------------
# Variant detection pipeline to output ALL raw SNPs.
# (LOCAL SCRIPT)
#***************************************************************************

# The user must specify the path to the repository, GATK and SnpEff
REPOSITORY_PATH="/path/to/repository"
GATK_PATH="/path/to/gatk"
SNPEFF_PATH="/path/to/snpeff"

python ./local/4_SelectFilterAnnotateVariants.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -snpeff-path $SNPEFF_PATH -population ALL -version Tcas3.30 -filter-suffix raw_SNP

