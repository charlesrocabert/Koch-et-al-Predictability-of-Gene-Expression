#!/usr/bin/env Rscript

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# 2_SplitBamMap.R
# ---------------
# Split the bam map into two files containing each complete version.
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

rm(list=ls())


##################
#      MAIN      #
##################

#--------------------------------#
# 1) Read command line arguments #
#--------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
setwd(REPOSITORY_PATH)

#--------------------------------#
# 2) Clean the BAM map           #
#--------------------------------#
bmap = read.table("data/tribolium_bam/bam_map.csv", h=T, sep=";")
bmap = bmap[!is.na(bmap$fitness) & !is.na(bmap$run_date) & !is.na(bmap$run_index),]

#--------------------------------#
# 3) Build Tcas3.30 BAM map      #
#--------------------------------#
bmap_tcas3_30 = bmap[bmap$annotation=="Version-2016-02-11",]
bmap_tcas3_30 = bmap_tcas3_30[!duplicated(bmap_tcas3_30$sample),]
write.table(bmap_tcas3_30, file="data/tribolium_bam/bam_map_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")

#--------------------------------#
# 4) Build Tcas5.2 BAM map       #
#--------------------------------#
bmap_tcas5_2 = bmap[bmap$annotation=="Version-2017-03-28",]
bmap_tcas5_2 = bmap_tcas5_2[!duplicated(bmap_tcas5_2$sample),]
write.table(bmap_tcas5_2, file="data/tribolium_bam/bam_map_Tcas5.2.csv", quote=F, row.names=F, col.names=T, sep=";")

#--------------------------------#
# 5) Save list of samples        #
#--------------------------------#
ALL   = bmap_tcas3_30[, seq(1,13)]
CT    = bmap_tcas3_30[bmap_tcas3_30$source_env=="CT", seq(1,13)]
HD    = bmap_tcas3_30[bmap_tcas3_30$source_env=="HD", seq(1,13)]
CT_HD = bmap_tcas3_30[bmap_tcas3_30$source_env=="CT" | bmap_tcas3_30$source_env=="HD", seq(1,13)]
LINES = c("L1", "L2", "L3", "L4", "L5", "L6", "Mx1", "Mx2", "Mx3")

### ALL ###
write.table(ALL, file="data/tribolium_bam/samples_ALL_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")
write.table(ALL[ALL$generation==1,], file="data/tribolium_bam/samples_ALL_G1_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")
write.table(ALL[ALL$generation==21,], file="data/tribolium_bam/samples_ALL_G21_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")

### CT ###
write.table(CT, file="data/tribolium_bam/samples_CT_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")
write.table(CT[CT$generation==1 & CT$line%in%LINES,], file="data/tribolium_bam/samples_CT_G1_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")
write.table(CT[CT$generation==21 & CT$line%in%LINES,], file="data/tribolium_bam/samples_CT_G21_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")
for(line in LINES)
{
  write.table(CT[CT$generation==1 & CT$line==line,], file=paste0("data/tribolium_bam/samples_CT_G1_",line,"_Tcas3.30.csv"), quote=F, row.names=F, col.names=T, sep=";")
  write.table(CT[CT$generation==21 & CT$line==line,], file=paste0("data/tribolium_bam/samples_CT_G21_",line,"_Tcas3.30.csv"), quote=F, row.names=F, col.names=T, sep=";")
}

### HD ###
write.table(HD, file="data/tribolium_bam/samples_HD_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")
write.table(HD[HD$generation==1 & HD$line%in%LINES,], file="data/tribolium_bam/samples_HD_G1_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")
write.table(HD[HD$generation==21 & HD$line%in%LINES,], file="data/tribolium_bam/samples_HD_G21_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")
for(line in LINES)
{
  write.table(HD[HD$generation==1 & HD$line==line,], file=paste0("data/tribolium_bam/samples_HD_G1_",line,"_Tcas3.30.csv"), quote=F, row.names=F, col.names=T, sep=";")
  write.table(HD[HD$generation==21 & HD$line==line,], file=paste0("data/tribolium_bam/samples_HD_G21_",line,"_Tcas3.30.csv"), quote=F, row.names=F, col.names=T, sep=";")
}

### CT/HD ###
write.table(CT_HD[CT_HD$line%in%LINES,], file="data/tribolium_bam/samples_CT_HD_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")
write.table(CT_HD[CT_HD$generation==1 & CT_HD$line%in%LINES,], file="data/tribolium_bam/samples_CT_HD_G1_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")
write.table(CT_HD[CT_HD$generation==21 & CT_HD$line%in%LINES,], file="data/tribolium_bam/samples_CT_HD_G21_Tcas3.30.csv", quote=F, row.names=F, col.names=T, sep=";")

