#!/bin/bash
# coding: utf-8

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# 1_CopyExpressionFiles.sh
# ------------------------
# Copy prepared expression files to the phenotypes folder.
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

VERSION="Tcas3.30"

cp ../../data/tribolium_counts/Tribolium_castaneum_CT_G1_Tcas3.30_read_counts_READY.txt ../../data/tribolium_phenotypes/Tribolium_castaneum_CT_G1_Tcas3.30_EXPRESSION.txt
gzip -c ../../data/tribolium_phenotypes/Tribolium_castaneum_CT_G1_Tcas3.30_EXPRESSION.txt > ../../data/tribolium_phenotypes/Tribolium_castaneum_CT_G1_Tcas3.30_expression.txt.gz

cp ../../data/tribolium_counts/Tribolium_castaneum_HD_G1_Tcas3.30_read_counts_READY.txt ../../data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_Tcas3.30_EXPRESSION.txt
gzip -c ../../data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_Tcas3.30_EXPRESSION.txt > ../../data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_Tcas3.30_EXPRESSION.txt.gz

