#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 2_ComputeAFCs.R
# ---------------
# Compute AFCs for a given environment.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")

### Build the AFC dataset ###
build_AFC_dataset <- function( AF, environment, keep_fixed_alleles )
{
  MERGED_DATA = c()
  LINES       = c()
  if (environment == "CT")
  {
    LINES = c("L1", "L3", "L5", "L6", "Mx1", "Mx3")
  } else if (environment == "HD")
  {
    LINES = c("L1", "L2", "L3", "L5", "L6", "Mx1", "Mx2")
  }
  for(line in LINES)
  {
    G1          = filter(AF, GENERATION==1 & LINE==line)[,c("ID", "POS", "CHROM", "LINE", "AF")]
    G21_AF      = filter(AF, GENERATION==21 & LINE==line)[,"AF"]
    D           = data.frame(G1, G21_AF)
    names(D)    = c("ID", "POS", "CHROM", "LINE", "AF_G1", "AF_G21")
    D$AFC       = D$AF_G21-D$AF_G1
    MERGED_DATA = rbind(MERGED_DATA, D)
  }
  if (!keep_fixed_alleles)
  {
    MERGED_DATA = filter(MERGED_DATA, AF_G1>0 & AF_G1<1 & AF_G21>0 & AF_G21<1)
  }
  MERGED_DATA = MERGED_DATA[order(MERGED_DATA$CHROM, MERGED_DATA$POS),]
  return(MERGED_DATA)
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
# 2) Compute AFCs                #
#--------------------------------#
KEEP_FIXED_ALLELES = TRUE
for (env in c("CT", "HD"))
{
  ### 2.1) Load allelic frequencies data ###
  AF = readRDS(paste0("./data/tribolium_afc/AF_",env,".rds"))
  
  ### 2.2) Build the AFC dataset ###
  AFC = build_AFC_dataset(AF, env, KEEP_FIXED_ALLELES)
  
  ### 2.3) Polarize SNPs to get MAFs at G1 (to fit Gemma output) ###
  POS              = which(AFC$AF_G1 > 0.5)
  AFC[POS,]$AFC    = -AFC[POS,]$AFC
  AFC[POS,]$AF_G1  = 1-AFC[POS,]$AF_G1
  AFC[POS,]$AF_G21 = 1-AFC[POS,]$AF_G21
  
  ### 2.4) Save the dataset ###
  saveRDS(AFC, file=paste0("./data/tribolium_afc/AFC_",env,".rds"))
  write.table(AFC, file=paste0("./data/tribolium_afc/AFC_",env,".csv"), sep=";", row.names=F, col.names=T, quote=F)
}

