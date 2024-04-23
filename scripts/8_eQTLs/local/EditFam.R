#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# EditFam.R
# ---------
# Edit .fam file to add full-sib families and phenotypes.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())


##################
#      MAIN      #
##################

#--------------------------------------#
# 1) Read command line arguments       #
#--------------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<5)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
POPULATION      = args[2]
VERSION         = args[3]
IMPUTATION      = args[4]
PHENOTYPE       = args[5]
setwd(REPOSITORY_PATH)

#--------------------------------------#
# 2) Loading samples                   #
#--------------------------------------#
samples                = read.table(paste0("./data/tribolium_bam/samples_",POPULATION,"_",VERSION,".csv"), h=T, sep=";")
samples$fem            = as.factor(samples$fem)
samples$line           = as.factor(samples$line)
samples$fullsib_family = as.factor(samples$fullsib_family)
samples$halfsib_family = as.factor(samples$halfsib_family)
samples$source_env     = as.factor(samples$source_env)
samples$target_env     = as.factor(samples$target_env)
samples$run_index      = as.factor(samples$run_index)
samples$batch_index    = as.factor(samples$batch_index)

#--------------------------------------#
# 3) Add family structure in .fam file #
#--------------------------------------#
FAM           = read.table(paste0("./data/tribolium_eqtl/gemma/",POPULATION,"_",VERSION,IMPUTATION,"_",PHENOTYPE,".fam"), sep="\t", check.names=F)
rownames(FAM) = FAM[,2]
FAM           = FAM[samples$sample,]
FAM[,1]       = samples$fem
FAM[,3]       = samples$halfsib_family
FAM[,4]       = samples$fullsib_family
FAM           = FAM[,1:5]

#--------------------------------------#
# 4) Loading phenotypes                #
#--------------------------------------#
if (PHENOTYPE == "EXPRESSION")
{
  PHENO      = read.table(paste0("./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_EXPRESSION.txt"), sep="\t", h=T, check.names=F)
  gene_id    = PHENO[,1]
  tPHENO     = t(PHENO[,2:dim(PHENO)[2]])
  tPHENO     = tPHENO[samples$sample,]
  FAM        = cbind(FAM, tPHENO)
  names(FAM) = c("FID", "IID", "Father", "Mother", "Sex", gene_id)
}
if (PHENOTYPE == "PLASTICITY")
{
  PHENO      = read.table(paste0("./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_PLASTICITY.txt"), sep="\t", h=T, check.names=F)
  gene_id    = PHENO[,1]
  tPHENO     = t(PHENO[,2:dim(PHENO)[2]])
  tPHENO     = tPHENO[samples$sample,]
  FAM        = cbind(FAM, tPHENO)
  names(FAM) = c("FID", "IID", "Father", "Mother", "Sex", gene_id)
}
if (PHENOTYPE == "NOISE")
{
  PHENO      = read.table(paste0("./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_NOISE.txt"), sep="\t", h=T, check.names=F)
  gene_id    = PHENO[,1]
  tPHENO     = t(PHENO[,2:dim(PHENO)[2]])
  tPHENO     = tPHENO[samples$sample,]
  FAM        = cbind(FAM, tPHENO)
  names(FAM) = c("FID", "IID", "Father", "Mother", "Sex", gene_id)
}
if (PHENOTYPE == "FITNESS")
{
  fitness = read.table(paste0("./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_FITNESS.txt"), sep="\t", h=T, check.names=F)
  fitness = fitness[samples$sample,"fitness"]
  FAM     = cbind(FAM, fitness)
  names(FAM) = c("FID", "IID", "Father", "Mother", "Sex", "fitness")
}
phenotypes = colnames(FAM)[6:dim(FAM)[2]]

#--------------------------------------#
# 5) Save edited .fam file             #
#--------------------------------------#
write.table(FAM, file=paste0("./data/tribolium_eqtl/gemma/",POPULATION,"_",VERSION,IMPUTATION,"_",PHENOTYPE,".fam"), quote=F, sep="\t", row.names=F, col.names=F)
write.table(phenotypes, file=paste0("./data/tribolium_eqtl/gemma/",POPULATION,"_",VERSION,IMPUTATION,"_",PHENOTYPE,".pheno"), quote=F, sep="\t", row.names=F, col.names=F)

