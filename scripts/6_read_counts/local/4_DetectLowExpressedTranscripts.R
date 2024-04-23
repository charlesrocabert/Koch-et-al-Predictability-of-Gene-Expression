#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 5_DetectLowExpressedTranscripts.R
# ---------------------------------
# Detect low expressed transcripts and save the list of expressed ones.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("edgeR")


##################
#      MAIN      #
##################

#-------------------------------------------#
# 1) Read command line arguments            #
#-------------------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
POPULATION      = args[2]
VERSION         = args[3]
SUFFIX          = args[4]
setwd(REPOSITORY_PATH)

#-------------------------------------------#
# 2) Loading samples                        #
#-------------------------------------------#
samples                = read.table(paste0("./data/tribolium_bam/samples_",POPULATION,"_",VERSION,".csv"), h=T, sep=";")
samples$fem            = as.factor(samples$fem)
samples$line           = as.factor(samples$line)
samples$fullsib_family = as.factor(samples$fullsib_family)
samples$halfsib_family = as.factor(samples$halfsib_family)
samples$source_env     = as.factor(samples$source_env)
samples$target_env     = as.factor(samples$target_env)
samples$run_index      = as.factor(samples$run_index)
samples$batch_index    = as.factor(samples$batch_index)

#-------------------------------------------#
# 3) Loading read counts                    #
#-------------------------------------------#
d           = read.table(paste0("./data/tribolium_counts/Tribolium_castaneum_",POPULATION,"_",VERSION,"_",SUFFIX,".txt"), sep="\t", h=T, check.names=F)
M           = as.matrix(d[,2:dim(d)[2]])
rownames(M) = d[,1]
M           = M[,samples$sample]

#-------------------------------------------#
# 4) Creating DGEList object                #
#-------------------------------------------#
dge = DGEList(M, remove.zeros=T)

#-------------------------------------------#
# 5) Filtering low read counts              #
#-------------------------------------------#
keep = rowSums(cpm(dge)>1) >= 2
dgek = dge[keep, , keep.lib.sizes=FALSE]

#-------------------------------------------#
# 6) Save the list of expressed transcripts #
#-------------------------------------------#
write.table(rownames(dgek$counts), file=paste0("./data/tribolium_counts/expressed_transcripts_",POPULATION,"_",VERSION,".txt"), quote=F, sep="\t", row.names=F, col.names=F)

