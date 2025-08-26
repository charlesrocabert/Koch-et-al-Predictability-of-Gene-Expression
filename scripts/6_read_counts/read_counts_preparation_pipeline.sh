#!/bin/bash
# coding: utf-8

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# read_counts_preparation_pipeline.sh
# -----------------------------------
# Read counts preparation pipeline.
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
VERSION="Tcas3.30"

### 1) Select populations from ALL read counts ###
python ../3_select_population/local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix read_counts -out-population CT_HD_G1 -out-suffix read_counts -version $VERSION -task COUNTS
python ../3_select_population/local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix read_counts -out-population CT_G1 -out-suffix read_counts -version $VERSION -task COUNTS
python ../3_select_population/local/1_SelectPopulation.py -repository-path $REPOSITORY_PATH -in-population ALL -in-suffix read_counts -out-population HD_G1 -out-suffix read_counts -version $VERSION -task COUNTS

### 2) Transform CT-HD-G1 read counts ###
Rscript ./local/3_TransformReadCounts.R $REPOSITORY_PATH CT_HD_G1 read_counts CT_HD_G1 read_counts_READY $VERSION

### 3) Select sub-populations from CT/HD G1 transformed data ###
python ../3_select_population/local/1_SelectSubPopulation.py -repository-path $REPOSITORY_PATH -in-population CT_HD_G1 -in-suffix read_counts_READY -out-population CT_G1 -out-suffix read_counts_transformed -version $VERSION -task COUNTS
python ../3_select_population/local/1_SelectSubPopulation.py -repository-path $REPOSITORY_PATH -in-population CT_HD_G1 -in-suffix read_counts_READY -out-population HD_G1 -out-suffix read_counts_transformed -version $VERSION -task COUNTS

### 4) Detect low expressed transcripts ###
Rscript ./local/4_DetectLowExpressedTranscripts.R $REPOSITORY_PATH CT_G1 read_counts $VERSION
Rscript ./local/4_DetectLowExpressedTranscripts.R $REPOSITORY_PATH HD_G1 read_counts $VERSION

### 5) Standardize read counts ###
Rscript ./local/5_StandardizeReadCounts.R $REPOSITORY_PATH CT_G1 read_counts_transformed CT_G1 read_counts_READY $VERSION
Rscript ./local/5_StandardizeReadCounts.R $REPOSITORY_PATH HD_G1 read_counts_transformed HD_G1 read_counts_READY $VERSION

