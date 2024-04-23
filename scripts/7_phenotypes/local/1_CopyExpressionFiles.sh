#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_CopyExpressionFiles.sh
# ------------------------
# Copy prepared expression files to the phenotypes folder.
# (LOCAL SCRIPT)
#***************************************************************************

VERSION="Tcas3.30"

cp ../../data/tribolium_counts/Tribolium_castaneum_CT_G1_Tcas3.30_read_counts_READY.txt ../../data/tribolium_phenotypes/Tribolium_castaneum_CT_G1_Tcas3.30_EXPRESSION.txt
gzip -c ../../data/tribolium_phenotypes/Tribolium_castaneum_CT_G1_Tcas3.30_EXPRESSION.txt > ../../data/tribolium_phenotypes/Tribolium_castaneum_CT_G1_Tcas3.30_expression.txt.gz

cp ../../data/tribolium_counts/Tribolium_castaneum_HD_G1_Tcas3.30_read_counts_READY.txt ../../data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_Tcas3.30_EXPRESSION.txt
gzip -c ../../data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_Tcas3.30_EXPRESSION.txt > ../../data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_Tcas3.30_EXPRESSION.txt.gz

