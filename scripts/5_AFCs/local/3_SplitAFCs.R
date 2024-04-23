#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 3_SplitAFCs.R
# -------------
# Split AFCs by line and environment.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")

### Pivot CT lines AFC values ###
pivot_AFC_CT_table <- function( AFC_CT, transform_AFCs )
{
  #------------------------------------------#
  # 1) Load L1 data                          #
  #------------------------------------------#
  L1_G1_HD     = filter(AFC_CT, LINE == "L1")[,c("ID","AF_G1")]
  L1_G21_HD    = filter(AFC_CT, LINE == "L1")[,c("ID","AF_G21")]
  L1_HD        = merge(L1_G1_HD, L1_G21_HD, by="ID")
  L1_HD$L1_AFC = L1_HD$AF_G21-L1_HD$AF_G1
  names(L1_HD) = c("ID", "L1_HD_G1", "L1_HD_G21", "L1_AFC")
  #------------------------------------------#
  # 2) Load L3 data                          #
  #------------------------------------------#
  L3_G1_HD     = filter(AFC_CT, LINE == "L3")[,c("ID","AF_G1")]
  L3_G21_HD    = filter(AFC_CT, LINE == "L3")[,c("ID","AF_G21")]
  L3_HD        = merge(L3_G1_HD, L3_G21_HD, by="ID")
  L3_HD$L3_AFC = L3_HD$AF_G21-L3_HD$AF_G1
  names(L3_HD) = c("ID", "L3_HD_G1", "L3_HD_G21", "L3_AFC")
  #------------------------------------------#
  # 3) Load L4 data                          #
  #------------------------------------------#
  L5_G1_HD     = filter(AFC_CT, LINE == "L5")[,c("ID","AF_G1")]
  L5_G21_HD    = filter(AFC_CT, LINE == "L5")[,c("ID","AF_G21")]
  L5_HD        = merge(L5_G1_HD, L5_G21_HD, by="ID")
  L5_HD$L5_AFC = L5_HD$AF_G21-L5_HD$AF_G1
  names(L5_HD) = c("ID", "L5_HD_G1", "L5_HD_G21", "L5_AFC")
  #------------------------------------------#
  # 4) Load L5 data                          #
  #------------------------------------------#
  L6_G1_HD     = filter(AFC_CT, LINE == "L6")[,c("ID","AF_G1")]
  L6_G21_HD    = filter(AFC_CT, LINE == "L6")[,c("ID","AF_G21")]
  L6_HD        = merge(L6_G1_HD, L6_G21_HD, by="ID")
  L6_HD$L6_AFC = L6_HD$AF_G21-L6_HD$AF_G1
  names(L6_HD) = c("ID", "L6_HD_G1", "L6_HD_G21", "L6_AFC")
  #-------------------------------------------#
  # 5) Load Mx1 data                          #
  #-------------------------------------------#
  Mx1_G1_HD      = filter(AFC_CT, LINE == "Mx1")[,c("ID","AF_G1")]
  Mx1_G21_HD     = filter(AFC_CT, LINE == "Mx1")[,c("ID","AF_G21")]
  Mx1_HD         = merge(Mx1_G1_HD, Mx1_G21_HD, by="ID")
  Mx1_HD$Mx1_AFC = Mx1_HD$AF_G21-Mx1_HD$AF_G1
  names(Mx1_HD) = c("ID", "Mx1_HD_G1", "Mx1_HD_G21", "Mx1_AFC")
  #-------------------------------------------#
  # 6) Load Mx3 data                          #
  #-------------------------------------------#
  Mx3_G1_HD      = filter(AFC_CT, LINE == "Mx3")[,c("ID","AF_G1")]
  Mx3_G21_HD     = filter(AFC_CT, LINE == "Mx3")[,c("ID","AF_G21")]
  Mx3_HD         = merge(Mx3_G1_HD, Mx3_G21_HD, by="ID")
  Mx3_HD$Mx3_AFC = Mx3_HD$AF_G21-Mx3_HD$AF_G1
  names(Mx3_HD) = c("ID", "Mx3_HD_G1", "Mx3_HD_G21", "Mx3_AFC")
  #-------------------------------------------#
  # 7) Build the AFC table                    #
  #-------------------------------------------#
  AFC_TABLE = merge(L1_HD, L3_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, L5_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, L6_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, Mx1_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, Mx3_HD, by="ID")
  #-------------------------------------------#
  # 8) Normalize AFs                          #
  #-------------------------------------------#
  coln = c("L1_HD_G1", "L1_HD_G21", "L3_HD_G1", "L3_HD_G21", "L5_HD_G1", "L5_HD_G21", "L6_HD_G1", "L6_HD_G21", "Mx1_HD_G1", "Mx1_HD_G21", "Mx3_HD_G1", "Mx3_HD_G21")
  if (transform_AFCs)
  {
    AFC_TABLE[,coln] = asin(sqrt(AFC_TABLE[,coln]))
    AFC_TABLE[,coln] = scale(AFC_TABLE[,coln], center=T, scale=T)
  }
  #-------------------------------------------#
  # 9) Split ID in pos and chr                #
  #-------------------------------------------#
  out  = strsplit(AFC_TABLE$ID,'-') 
  out2 = do.call(rbind, out)
  AFC_TABLE$CHR = out2[,1]
  AFC_TABLE$POS = out2[,2]
  #-------------------------------------------#
  # 10) Save absolute values of AFCs per line #
  #-------------------------------------------#
  AFC_TABLE$L1_AFC_abs  = abs(AFC_TABLE$L1_AFC)
  AFC_TABLE$L3_AFC_abs  = abs(AFC_TABLE$L3_AFC)
  AFC_TABLE$L5_AFC_abs  = abs(AFC_TABLE$L5_AFC)
  AFC_TABLE$L6_AFC_abs  = abs(AFC_TABLE$L6_AFC)
  AFC_TABLE$Mx1_AFC_abs = abs(AFC_TABLE$Mx1_AFC)
  AFC_TABLE$Mx3_AFC_abs = abs(AFC_TABLE$Mx3_AFC)
  #-------------------------------------------#
  # 12) Return the table                      #
  #-------------------------------------------#
  return(AFC_TABLE)
}

### Pivot HD lines AFC values ###
pivot_AFC_HD_table <- function( AFC_HD, transform_AFCs )
{
  #------------------------------------------#
  # 1) Load L1 data                          #
  #------------------------------------------#
  L1_G1_HD     = filter(AFC_HD, LINE == "L1")[,c("ID","AF_G1")]
  L1_G21_HD    = filter(AFC_HD, LINE == "L1")[,c("ID","AF_G21")]
  L1_HD        = merge(L1_G1_HD, L1_G21_HD, by="ID")
  L1_HD$L1_AFC = L1_HD$AF_G21-L1_HD$AF_G1
  names(L1_HD) = c("ID", "L1_HD_G1", "L1_HD_G21", "L1_AFC")
  #------------------------------------------#
  # 2) Load L2 data                          #
  #------------------------------------------#
  L2_G1_HD     = filter(AFC_HD, LINE == "L2")[,c("ID","AF_G1")]
  L2_G21_HD    = filter(AFC_HD, LINE == "L2")[,c("ID","AF_G21")]
  L2_HD        = merge(L2_G1_HD, L2_G21_HD, by="ID")
  L2_HD$L2_AFC = L2_HD$AF_G21-L2_HD$AF_G1
  names(L2_HD) = c("ID", "L2_HD_G1", "L2_HD_G21", "L2_AFC")
  #------------------------------------------#
  # 3) Load L3 data                          #
  #------------------------------------------#
  L3_G1_HD     = filter(AFC_HD, LINE == "L3")[,c("ID","AF_G1")]
  L3_G21_HD    = filter(AFC_HD, LINE == "L3")[,c("ID","AF_G21")]
  L3_HD        = merge(L3_G1_HD, L3_G21_HD, by="ID")
  L3_HD$L3_AFC = L3_HD$AF_G21-L3_HD$AF_G1
  names(L3_HD) = c("ID", "L3_HD_G1", "L3_HD_G21", "L3_AFC")
  #------------------------------------------#
  # 4) Load L4 data                          #
  #------------------------------------------#
  L5_G1_HD     = filter(AFC_HD, LINE == "L5")[,c("ID","AF_G1")]
  L5_G21_HD    = filter(AFC_HD, LINE == "L5")[,c("ID","AF_G21")]
  L5_HD        = merge(L5_G1_HD, L5_G21_HD, by="ID")
  L5_HD$L5_AFC = L5_HD$AF_G21-L5_HD$AF_G1
  names(L5_HD) = c("ID", "L5_HD_G1", "L5_HD_G21", "L5_AFC")
  #------------------------------------------#
  # 5) Load L5 data                          #
  #------------------------------------------#
  L6_G1_HD     = filter(AFC_HD, LINE == "L6")[,c("ID","AF_G1")]
  L6_G21_HD    = filter(AFC_HD, LINE == "L6")[,c("ID","AF_G21")]
  L6_HD        = merge(L6_G1_HD, L6_G21_HD, by="ID")
  L6_HD$L6_AFC = L6_HD$AF_G21-L6_HD$AF_G1
  names(L6_HD) = c("ID", "L6_HD_G1", "L6_HD_G21", "L6_AFC")
  #-------------------------------------------#
  # 6) Load Mx1 data                          #
  #-------------------------------------------#
  Mx1_G1_HD      = filter(AFC_HD, LINE == "Mx1")[,c("ID","AF_G1")]
  Mx1_G21_HD     = filter(AFC_HD, LINE == "Mx1")[,c("ID","AF_G21")]
  Mx1_HD         = merge(Mx1_G1_HD, Mx1_G21_HD, by="ID")
  Mx1_HD$Mx1_AFC = Mx1_HD$AF_G21-Mx1_HD$AF_G1
  names(Mx1_HD) = c("ID", "Mx1_HD_G1", "Mx1_HD_G21", "Mx1_AFC")
  #-------------------------------------------#
  # 7) Load Mx2 data                          #
  #-------------------------------------------#
  Mx2_G1_HD      = filter(AFC_HD, LINE == "Mx2")[,c("ID","AF_G1")]
  Mx2_G21_HD     = filter(AFC_HD, LINE == "Mx2")[,c("ID","AF_G21")]
  Mx2_HD         = merge(Mx2_G1_HD, Mx2_G21_HD, by="ID")
  Mx2_HD$Mx2_AFC = Mx2_HD$AF_G21-Mx2_HD$AF_G1
  names(Mx2_HD) = c("ID", "Mx2_HD_G1", "Mx2_HD_G21", "Mx2_AFC")
  #-------------------------------------------#
  # 8) Build the AFC table                    #
  #-------------------------------------------#
  AFC_TABLE = merge(L1_HD, L2_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, L3_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, L5_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, L6_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, Mx1_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, Mx2_HD, by="ID")
  #-------------------------------------------#
  # 9) Normalize AFs                          #
  #-------------------------------------------#
  coln = c("L1_HD_G1", "L1_HD_G21", "L2_HD_G1", "L2_HD_G21", "L3_HD_G1", "L3_HD_G21", "L5_HD_G1", "L5_HD_G21", "L6_HD_G1", "L6_HD_G21", "Mx1_HD_G1", "Mx1_HD_G21", "Mx2_HD_G1", "Mx2_HD_G21")
  if (transform_AFCs)
  {
    AFC_TABLE[,coln] = asin(sqrt(AFC_TABLE[,coln]))
    AFC_TABLE[,coln] = scale(AFC_TABLE[,coln], center=T, scale=T)
  }
  #-------------------------------------------#
  # 10) Split ID in pos and chr               #
  #-------------------------------------------#
  out  = strsplit(AFC_TABLE$ID,'-') 
  out2 = do.call(rbind, out)
  AFC_TABLE$CHR = out2[,1]
  AFC_TABLE$POS = out2[,2]
  #-------------------------------------------#
  # 11) Save absolute values of AFCs per line #
  #-------------------------------------------#
  AFC_TABLE$L1_AFC_abs  = abs(AFC_TABLE$L1_AFC)
  AFC_TABLE$L2_AFC_abs  = abs(AFC_TABLE$L2_AFC)
  AFC_TABLE$L3_AFC_abs  = abs(AFC_TABLE$L3_AFC)
  AFC_TABLE$L5_AFC_abs  = abs(AFC_TABLE$L5_AFC)
  AFC_TABLE$L6_AFC_abs  = abs(AFC_TABLE$L6_AFC)
  AFC_TABLE$Mx1_AFC_abs = abs(AFC_TABLE$Mx1_AFC)
  AFC_TABLE$Mx2_AFC_abs = abs(AFC_TABLE$Mx2_AFC)
  #-------------------------------------------#
  # 12) Return the table                      #
  #-------------------------------------------#
  return(AFC_TABLE)
}

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
# 2) Load the AFC values         #
#--------------------------------#
AFC_CT       = readRDS("./data/tribolium_afc/AFC_CT.rds")
AFC_HD       = readRDS("./data/tribolium_afc/AFC_HD.rds")
AFC_CT_TABLE = pivot_AFC_CT_table(AFC_CT, FALSE)
AFC_HD_TABLE = pivot_AFC_HD_table(AFC_HD, FALSE)

#--------------------------------#
# 3) Save the dataset            #
#--------------------------------#
saveRDS(AFC_CT_TABLE$L1_AFC_abs, "./data/tribolium_afc/vectors/AFC_CT_L1.rds")
saveRDS(AFC_CT_TABLE$L3_AFC_abs, "./data/tribolium_afc/vectors/AFC_CT_L3.rds")
saveRDS(AFC_CT_TABLE$L5_AFC_abs, "./data/tribolium_afc/vectors/AFC_CT_L5.rds")
saveRDS(AFC_CT_TABLE$L6_AFC_abs, "./data/tribolium_afc/vectors/AFC_CT_L6.rds")
saveRDS(AFC_CT_TABLE$Mx1_AFC_abs, "./data/tribolium_afc/vectors/AFC_CT_Mx1.rds")
saveRDS(AFC_CT_TABLE$Mx3_AFC_abs, "./data/tribolium_afc/vectors/AFC_CT_Mx3.rds")

saveRDS(AFC_HD_TABLE$L1_AFC_abs, "./data/tribolium_afc/vectors/AFC_HD_L1.rds")
saveRDS(AFC_HD_TABLE$L2_AFC_abs, "./data/tribolium_afc/vectors/AFC_HD_L2.rds")
saveRDS(AFC_HD_TABLE$L3_AFC_abs, "./data/tribolium_afc/vectors/AFC_HD_L3.rds")
saveRDS(AFC_HD_TABLE$L5_AFC_abs, "./data/tribolium_afc/vectors/AFC_HD_L5.rds")
saveRDS(AFC_HD_TABLE$L6_AFC_abs, "./data/tribolium_afc/vectors/AFC_HD_L6.rds")
saveRDS(AFC_HD_TABLE$Mx1_AFC_abs, "./data/tribolium_afc/vectors/AFC_HD_Mx1.rds")
saveRDS(AFC_HD_TABLE$Mx2_AFC_abs, "./data/tribolium_afc/vectors/AFC_HD_Mx2.rds")


