#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 4_ComputeRelativeFitness.R
# --------------------------
# Compute relative fitness per line for CT and HD individuals (G1).
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
if (length(args)<2)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
POPULATION      = args[2]
VERSION         = args[3]
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
# 3) Compute relative fitness per line #
#--------------------------------------#
ENVIRONMENTS = c("CT", "HD")
SAMPLES_vec  = c()
FITNESS_vec  = c()
for (env in ENVIRONMENTS)
{
  SAMPLES     = samples[samples$target_env==env,]
  mean_w      = mean(SAMPLES$fitness)
  w           = SAMPLES$fitness/mean_w
  SAMPLES_vec = c(SAMPLES_vec, SAMPLES$sample)
  FITNESS_vec = c(FITNESS_vec, w)
}
FITNESS_DATA           = data.frame(SAMPLES_vec, FITNESS_vec)
names(FITNESS_DATA)    = c("sample", "fitness")
rownames(FITNESS_DATA) = FITNESS_DATA$sample

CT_FITNESS = FITNESS_DATA[samples[samples$target_env=="CT","sample"],]
HD_FITNESS = FITNESS_DATA[samples[samples$target_env=="HD","sample"],]

#--------------------------------------#
# 4) Correct for fullsib effect        #
#--------------------------------------#
# reg = lm(CT_FITNESS$fitness~samples[samples$target_env=="CT","fem"])
# CT_FITNESS$fitness_adjusted = reg$residuals
#
# reg = lm(HD_FITNESS$fitness~samples[samples$target_env=="HD","fem"])
# HD_FITNESS$fitness_adjusted = reg$residuals

#--------------------------------------#
# 5) Quantile normalization            #
#--------------------------------------#
mat                = CT_FITNESS$fitness
mat                = rank(mat, ties.method="average")
mat                = qnorm(mat/(nrow(CT_FITNESS)+1))
CT_FITNESS$fitness = mat
rm(mat)

# mat                         = CT_FITNESS$fitness_adjusted
# mat                         = rank(mat, ties.method="average")
# mat                         = qnorm(mat/(nrow(CT_FITNESS)+1))
# CT_FITNESS$fitness_adjusted = mat
# rm(mat)

mat                = HD_FITNESS$fitness
mat                = rank(mat, ties.method="average")
mat                = qnorm(mat/(nrow(HD_FITNESS)+1))
HD_FITNESS$fitness = mat
rm(mat)

# mat                         = HD_FITNESS$fitness_adjusted
# mat                         = rank(mat, ties.method="average")
# mat                         = qnorm(mat/(nrow(HD_FITNESS)+1))
# HD_FITNESS$fitness_adjusted = mat
# rm(mat)

#--------------------------------------#
# 5) Save the resulting datasets       #
#--------------------------------------#
write.table(CT_FITNESS, file=paste0("./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_FITNESS.txt"), quote=F, sep="\t", row.names=T, col.names=T)
system(paste0("gzip -c ./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_FITNESS.txt > ./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_FITNESS.txt.gz"))

write.table(HD_FITNESS, file=paste0("./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_FITNESS.txt"), quote=F, sep="\t", row.names=T, col.names=T)
system(paste0("gzip -c ./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_FITNESS.txt > ./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_FITNESS.txt.gz"))

