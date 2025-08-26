#!/bin/bash
# coding: utf-8

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# 6_merge_eQTLs_pipeline.sh
# -------------------------
# Merge eQTL associations.
# (LOCAL SCRIPT)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*******************************************************************************


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

