#!/bin/bash
# coding: utf-8

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# variant_detection_pipeline.sh
# -----------------------------
# Variant detection pipeline to output ALL raw SNPs.
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

# The user must specify the path to the repository, GATK and SnpEff
REPOSITORY_PATH="/path/to/repository"
GATK_PATH="/path/to/gatk"
SNPEFF_PATH="/path/to/snpeff"

python ./local/4_SelectFilterAnnotateVariants.py -repository-path $REPOSITORY_PATH -gatk-path $GATK_PATH -snpeff-path $SNPEFF_PATH -population ALL -version Tcas3.30 -filter-suffix raw_SNP

