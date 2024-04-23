#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 6_StandardizeReadCounts.R
# -------------------------
# Standardize a gene expression dataset by quantile normalization.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("edgeR")
library("limma")
library("ggfortify")


##################
#      MAIN      #
##################

#---------------------------------------------#
# 1) Read command line arguments              #
#---------------------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<6)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
IN_POPULATION   = args[2]
IN_SUFFIX       = args[3]
OUT_POPULATION  = args[4]
OUT_SUFFIX      = args[5]
VERSION         = args[6]

stopifnot(paste0(IN_POPULATION,"_",IN_SUFFIX) != paste0(OUT_POPULATION,"_",OUT_SUFFIX))
setwd(REPOSITORY_PATH)

#---------------------------------------------#
# 2) Loading samples                          #
#---------------------------------------------#
samples                = read.table(paste0("./data/tribolium_bam/samples_",OUT_POPULATION,"_",VERSION,".csv"), h=T, sep=";")
samples$fem            = as.factor(samples$fem)
samples$line           = as.factor(samples$line)
samples$fullsib_family = as.factor(samples$fullsib_family)
samples$halfsib_family = as.factor(samples$halfsib_family)
samples$source_env     = as.factor(samples$source_env)
samples$target_env     = as.factor(samples$target_env)
samples$run_index      = as.factor(samples$run_index)
samples$batch_index    = as.factor(samples$batch_index)

#---------------------------------------------#
# 3) Loading list of expressed transcripts    #
#---------------------------------------------#
transcripts = read.table(paste0("./data/tribolium_counts/expressed_transcripts_",IN_POPULATION,"_",VERSION,".txt"), sep="\t", h=T, check.names=F)[,1]

#---------------------------------------------#
# 4) Loading read counts                      #
#---------------------------------------------#
d       = read.table(paste0("./data/tribolium_counts/Tribolium_castaneum_",IN_POPULATION,"_",VERSION,"_",IN_SUFFIX,".txt"), sep="\t", h=T, check.names=F)
d       = d[d$gene_id%in%transcripts,]
gene_id = d[,1]
X       = d[,samples$sample]

#---------------------------------------------#
# 5) Quantile normalization                   #
#---------------------------------------------#
for(i in 1:dim(X)[1])
{
  if(i%%1000==0)
  {
    print(i)
  }
  mat   = X[i,]
  mat   = rank(mat, ties.method="average")
  mat   = qnorm(mat/(ncol(X)+1))
  X[i,] = mat
}
rm(i, mat)
X = as.data.frame(X)
X = cbind(gene_id, X)

#---------------------------------------------#
# 6) Save the normalized and adjusted dataset #
#---------------------------------------------#
filename = paste0("./data/tribolium_counts/Tribolium_castaneum_",OUT_POPULATION,"_",VERSION,"_",OUT_SUFFIX,".txt")
write.table(X, file=filename, quote=F, sep="\t", row.names=F, col.names=T)
system(paste0("gzip -c ",filename," > ",filename,".gz"))

