#!/bin/bash
# coding: utf-8

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# VCF_HD_G1_eQTL_imputed_pipeline.sh
# ----------------------------------
# SNP pipeline for eQTL (HD G1 imputed population).
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
GATK_PATH="/path/to/gatk"

python ./local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix imputed -out-population HD_G1 -out-suffix imputed -version Tcas3.30 -task VCF
python ./local/2_FilterGenotypes.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population HD_G1 -version Tcas3.30 -imputed -suffix eQTL_imputed
python ./local/3_VariantsToTable.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -population HD_G1 -version Tcas3.30 -suffix eQTL_imputed -annotation

