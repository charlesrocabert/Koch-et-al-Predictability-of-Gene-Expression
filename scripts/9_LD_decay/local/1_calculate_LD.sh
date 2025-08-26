#!/bin/bash
# coding: utf-8

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# 1_calculate_LD.sh
# -----------------
# Calculate pairwise LD (r2) for each chromosome.
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


cd /Users/charlesrocabert/git/charlesrocabert/Tribolium-Polygenic-Adaptation/

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1) Define main parameters           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
VERSION="Tcas3.30"
SUFFIX="raw_SNP"
PLINK="/Users/charlesrocabert/plink1.9/plink"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2) Calculate LD for each chromosome #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Calculate LD for each chromosome"

for CHR in ChLG2 ChLG3 ChLG4 ChLG5 ChLG6 ChLG7 ChLG8 ChLG9 ChLGX
do
    echo "   > Processing chromosome $CHR"
    $PLINK --vcf ./data/tribolium_snp/Tribolium_castaneum_ALL_$VERSION\_$SUFFIX.vcf.gz \
    --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --maf 0.05 --geno 0.5 --mind 0.5 --chr $CHR \
    --thin 0.01 -r2 gz --ld-window 100 --ld-window-kb 1000 \
    --ld-window-r2 0 \
    --make-bed --out ./data/tribolium_ld/$CHR
    rm ./data/tribolium_ld/$CHR.nosex
    rm ./data/tribolium_ld/$CHR.log
    rm ./data/tribolium_ld/$CHR.bim
    rm ./data/tribolium_ld/$CHR.bed
    rm ./data/tribolium_ld/$CHR.fam
    rm ./data/tribolium_ld/$CHR.irem
    gzip -d ./data/tribolium_ld/$CHR.ld.gz
done

