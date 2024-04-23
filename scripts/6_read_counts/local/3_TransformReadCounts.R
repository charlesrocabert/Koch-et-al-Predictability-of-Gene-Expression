#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 4_TransformReadCounts.R
# -----------------------
# Transform a gene expression dataset given the user options.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("edgeR")
library("limma")
library("ggfortify")
library("tibble")


##################
#      MAIN      #
##################

#--------------------------------------------------#
# 1) Read command line arguments                   #
#--------------------------------------------------#
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

#--------------------------------------------------#
# 2) Loading samples                               #
#--------------------------------------------------#
samples                = read.table(paste0("./data/tribolium_bam/samples_",OUT_POPULATION,"_",VERSION,".csv"), h=T, sep=";")
samples$fem            = as.factor(samples$fem)
samples$line           = as.factor(samples$line)
samples$fullsib_family = as.factor(samples$fullsib_family)
samples$halfsib_family = as.factor(samples$halfsib_family)
samples$source_env     = as.factor(samples$source_env)
samples$target_env     = as.factor(samples$target_env)
samples$run_index      = as.factor(samples$run_index)
samples$batch_index    = as.factor(samples$batch_index)

#--------------------------------------------------#
# 3) Loading read counts                           #
#--------------------------------------------------#
d           = read.table(paste0("./data/tribolium_counts/Tribolium_castaneum_",IN_POPULATION,"_",VERSION,"_",IN_SUFFIX,".txt"), sep="\t", h=T, check.names=F)
M           = as.matrix(d[,2:dim(d)[2]])
rownames(M) = d[,1]
M           = M[,samples$sample]
X           = DGEList(M, remove.zeros=T, group=samples$target_env)

#--------------------------------------------------#
# 4) Filtering                                     #
#--------------------------------------------------#
keep = rowSums(cpm(X)>1) >=2
X    = X[keep, , keep.lib.sizes=FALSE]

#--------------------------------------------------#
# 5) TMM normalization                             #
#--------------------------------------------------#
X = calcNormFactors(X, method="TMM")

#--------------------------------------------------#
# 6) Remove sequencing run, batch and line effects #
#--------------------------------------------------#
design = model.matrix(~0 + batch_index+run_index+target_env, data=samples)
V1     = voom(X, design)
corfit = duplicateCorrelation(V1, design, block=as.integer(samples$line))
V2     = voom(X, design, block=samples$line, correlation=corfit$consensus)
X      = as.data.frame(V2$E)
X      = rownames_to_column(X, "gene_id")

#--------------------------------------------------#
# 7) Save the normalized and adjusted dataset      #
#--------------------------------------------------#
filename = paste0("./data/tribolium_counts/Tribolium_castaneum_",OUT_POPULATION,"_",VERSION,"_",OUT_SUFFIX,".txt")
write.table(X, file=filename, quote=F, sep="\t", row.names=F, col.names=T)
system(paste0("gzip -c ",filename," > ",filename,".gz"))

