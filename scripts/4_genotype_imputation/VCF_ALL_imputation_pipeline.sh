#!/bin/bash
# coding: utf-8

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# VCF_ALL_imputation_pipeline.sh
# ------------------------------
# Genotypes imputation pipeline (ALL population).
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

# The user must specify the path to the repository and Beagle
REPOSITORY_PATH="/path/to/repository"
BEAGLE_PATH="/path/to/beagle"

python ./local/3_ImputeGenotypes.py -repository-path $REPOSITORY_PATH -beagle $BEAGLE_PATH -population ALL -version Tcas3.30 -suffix Beagle

