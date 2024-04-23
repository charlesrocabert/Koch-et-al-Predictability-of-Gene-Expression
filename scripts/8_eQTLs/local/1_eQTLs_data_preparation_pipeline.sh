#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_eQTLs_data_preparation_pipeline.sh
# ------------------------------------
# eQTLs data preparation pipeline (before eQTLs mapping).
# (LOCAL SCRIPT)
#***************************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1) Define main parameters  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# The user must specify the path to the repository and GATK
REPOSITORY_PATH="/path/to/repository"
VERSION="Tcas3.30"
IMPUTATION="_imputed"
PLINK2="/Users/charlesrocabert/plink2/plink2"

cd $REPOSITORY_PATH

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2) Build gemma input files #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Generate .bed and .bim files"
$PLINK2 --vcf ./data/tribolium_snp/Tribolium_castaneum_CT_G1_$VERSION\_eQTL$IMPUTATION.vcf --allow-extra-chr --make-bed --out ./data/tribolium_eqtl/gemma/CT_G1_$VERSION$IMPUTATION\_EXPRESSION
$PLINK2 --vcf ./data/tribolium_snp/Tribolium_castaneum_CT_G1_$VERSION\_eQTL$IMPUTATION.vcf --allow-extra-chr --make-bed --out ./data/tribolium_eqtl/gemma/CT_G1_$VERSION$IMPUTATION\_FITNESS
$PLINK2 --vcf ./data/tribolium_snp/Tribolium_castaneum_HD_G1_$VERSION\_eQTL$IMPUTATION.vcf --allow-extra-chr --make-bed --out ./data/tribolium_eqtl/gemma/HD_G1_$VERSION$IMPUTATION\_EXPRESSION
$PLINK2 --vcf ./data/tribolium_snp/Tribolium_castaneum_HD_G1_$VERSION\_eQTL$IMPUTATION.vcf --allow-extra-chr --make-bed --out ./data/tribolium_eqtl/gemma/HD_G1_$VERSION$IMPUTATION\_FITNESS

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3) Edit Fam files          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Edit .fam files"
Rscript ./scripts/8_eQTLs/local/EditFam.R $REPOSITORY_PATH CT_G1 $VERSION $IMPUTATION EXPRESSION
Rscript ./scripts/8_eQTLs/local/EditFam.R $REPOSITORY_PATH CT_G1 $VERSION $IMPUTATION FITNESS
Rscript ./scripts/8_eQTLs/local/EditFam.R $REPOSITORY_PATH HD_G1 $VERSION $IMPUTATION EXPRESSION
Rscript ./scripts/8_eQTLs/local/EditFam.R $REPOSITORY_PATH HD_G1 $VERSION $IMPUTATION FITNESS

