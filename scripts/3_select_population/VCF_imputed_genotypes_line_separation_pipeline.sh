#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# VCF_imputed_genotypes_line_separation_pipeline.sh
# -------------------------------------------------
# Imputed genotypes extraction for each env/generation/line.
# (LOCAL SCRIPT)
#***************************************************************************

# The user must specify the path to the repository and GATK
REPOSITORY_PATH="/path/to/repository"
GATK_PATH="/path/to/gatk"
DATA_PATH="../../data/tribolium_snp"

#---------------------------------------#
# 1) Split populations at generation 1  #
#---------------------------------------#

### 1.1) Create pooled population at CT-G1 ###
python ./local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix imputed -out-population CT_G1 -out-suffix imputed -version Tcas3.30 -task VCF
python ./local/3_VariantsToTable.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population CT_G1 -version Tcas3.30 -suffix imputed -annotation
rm $DATA_PATH/Tribolium_castaneum_CT_G1_Tcas3.30_imputed.vcf
mv $DATA_PATH/Tribolium_castaneum_CT_G1_Tcas3.30_imputed.vcf.gz $DATA_PATH/imputed/Tribolium_castaneum_CT_G1_Tcas3.30_imputed.vcf.gz
mv $DATA_PATH/snp_table_CT_G1_Tcas3.30_imputed.csv $DATA_PATH/imputed/snp_table_CT_G1_Tcas3.30_imputed.csv

### 1.2) Split availables lines at CT-G1 ###
declare -a LIN=("L1" "L3" "L5" "L6")
for lin in ${LIN[@]}; do
  echo $lin
  python ./local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix imputed -out-population CT_G1_$lin -out-suffix imputed -version Tcas3.30 -task VCF
  python ./local/3_VariantsToTable.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population CT_G1_$lin -version Tcas3.30 -suffix imputed -annotation
  rm $DATA_PATH/Tribolium_castaneum_CT_G1_$lin\_Tcas3.30_imputed.vcf
  mv $DATA_PATH/Tribolium_castaneum_CT_G1_$lin\_Tcas3.30_imputed.vcf.gz $DATA_PATH/imputed/Tribolium_castaneum_CT_G1_$lin\_Tcas3.30_imputed.vcf.gz
  mv $DATA_PATH/snp_table_CT_G1_$lin\_Tcas3.30_imputed.csv $DATA_PATH/imputed/snp_table_CT_G1_$lin\_Tcas3.30_imputed.csv
done

### 1.3) Create pooled population at HD-G1 ###
python ./local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix imputed -out-population HD_G1 -out-suffix imputed -version Tcas3.30 -task VCF
python ./local/3_VariantsToTable.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population HD_G1 -version Tcas3.30 -suffix imputed -annotation
rm $DATA_PATH/Tribolium_castaneum_HD_G1_Tcas3.30_imputed.vcf
mv $DATA_PATH/Tribolium_castaneum_HD_G1_Tcas3.30_imputed.vcf.gz $DATA_PATH/imputed/Tribolium_castaneum_HD_G1_Tcas3.30_imputed.vcf.gz
mv $DATA_PATH/snp_table_HD_G1_Tcas3.30_imputed.csv $DATA_PATH/imputed/snp_table_HD_G1_Tcas3.30_imputed.csv

### 1.4) Split availables lines at HD-G1 ###
declare -a LIN=("L1" "L3" "L5" "L6")
for lin in ${LIN[@]}; do
  echo $lin
  python ./local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix imputed -out-population HD_G1_$lin -out-suffix imputed -version Tcas3.30 -task VCF
  python ./local/3_VariantsToTable.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population HD_G1_$lin -version Tcas3.30 -suffix imputed -annotation
  rm $DATA_PATH/Tribolium_castaneum_HD_G1_$lin\_Tcas3.30_imputed.vcf
  mv $DATA_PATH/Tribolium_castaneum_HD_G1_$lin\_Tcas3.30_imputed.vcf.gz $DATA_PATH/imputed/Tribolium_castaneum_HD_G1_$lin\_Tcas3.30_imputed.vcf.gz
  mv $DATA_PATH/snp_table_HD_G1_$lin\_Tcas3.30_imputed.csv $DATA_PATH/imputed/snp_table_HD_G1_$lin\_Tcas3.30_imputed.csv
done

#---------------------------------------#
# 2) Split populations at generation 21 #
#---------------------------------------#

### 2.1) Split availables lines at CT-G21 ###
declare -a LIN=("L1" "L2" "L3" "L5" "L6" "Mx1" "Mx3")
for lin in ${LIN[@]}; do
  echo $lin
  python ./local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix imputed -out-population CT_G21_$lin -out-suffix imputed -version Tcas3.30 -task VCF
  python ./local/3_VariantsToTable.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population CT_G21_$lin -version Tcas3.30 -suffix imputed -annotation
  rm $DATA_PATH/Tribolium_castaneum_CT_G21_$lin\_Tcas3.30_imputed.vcf
  mv $DATA_PATH/Tribolium_castaneum_CT_G21_$lin\_Tcas3.30_imputed.vcf.gz $DATA_PATH/imputed/Tribolium_castaneum_CT_G21_$lin\_Tcas3.30_imputed.vcf.gz
  mv $DATA_PATH/snp_table_CT_G21_$lin\_Tcas3.30_imputed.csv $DATA_PATH/imputed/snp_table_CT_G21_$lin\_Tcas3.30_imputed.csv
done

### 2.2) Split availables lines at HD-G21 ###
declare -a LIN=("L1" "L2" "L3" "L5" "L6" "Mx1" "Mx2")
for lin in ${LIN[@]}; do
  echo $lin
  python ./local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix imputed -out-population HD_G21_$lin -out-suffix imputed -version Tcas3.30 -task VCF
  python ./local/3_VariantsToTable.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population HD_G21_$lin -version Tcas3.30 -suffix imputed -annotation
  rm $DATA_PATH/Tribolium_castaneum_HD_G21_$lin\_Tcas3.30_imputed.vcf
  mv $DATA_PATH/Tribolium_castaneum_HD_G21_$lin\_Tcas3.30_imputed.vcf.gz $DATA_PATH/imputed/Tribolium_castaneum_HD_G21_$lin\_Tcas3.30_imputed.vcf.gz
  mv $DATA_PATH/snp_table_HD_G21_$lin\_Tcas3.30_imputed.csv $DATA_PATH/imputed/snp_table_HD_G21_$lin\_Tcas3.30_imputed.csv
done

