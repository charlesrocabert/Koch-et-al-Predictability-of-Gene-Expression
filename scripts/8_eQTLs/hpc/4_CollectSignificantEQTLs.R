#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 4_CollectSignificantEQTLs.R
# ---------------------------
# Collect all significant eQTL associations given a p-value threshold.
# - Genomic correction is applied and the FDR is calculated,
# - eQTLs with a p-value< 0.05 are selected.
# (HPC SCRIPT --> run wrapper)
#***************************************************************************

rm(list=ls())

### Upload a file to Allas ###
upload_to_allas <- function( filename, bucket )
{
  cmdline = paste0("rclone copyto ", filename, " allas:", bucket, "/eqtl/", filename)
  system(cmdline)
}

### Calculate p-value inflation lambda ###
calculate_lambda <- function( p )
{
  p      = p[!is.na(p)]
  n      = length(p)
  x2obs  = qchisq(p, 1, lower.tail=FALSE)
  x2exp  = qchisq(1:n/n, 1, lower.tail=FALSE)
  lambdA = median(x2obs)/median(x2exp)
  return(lambdA)
}

### Apply genomic control adjustment on p-values ###
calculate_gc <- function( p )
{
  n      = length(p)
  x2obs  = qchisq(p, 1, lower.tail=FALSE)
  x2exp  = qchisq(1:n/n, 1, lower.tail=FALSE)
  lambdA = median(x2obs)/median(x2exp)
  x2new  = x2obs/lambdA
  gc     = pchisq(x2new, df=1, lower.tail=FALSE)
  return(gc)
}

### Calculate the FDR ###
calculate_fdr <- function( p )
{
  fdr = p.adjust(p, method="fdr")
  return(fdr)
}


##################
#      MAIN      #
##################

#--------------------------------#
# 1) Read command line arguments #
#--------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
DATASET = args[1]
PATH    = args[2]

print(paste0(">>> Working with dataset ",DATASET,"..."))
setwd(PATH)

#--------------------------------#
# 2) Load the list of phenotypes #
#--------------------------------#
PHENOTYPES = read.table(paste0("/scratch/project_XXXXXXX/Tribolium-Polygenic-Adaptation/data/tribolium_eqtl/gemma/",DATASET,".pheno"), h=F, check.names=F)[,1]
N          = length(PHENOTYPES)

#--------------------------------#
# 3) Collect significant SNPs    #
#--------------------------------#
COLLECTED_EQTL      = c()
COLLECTED_PHENOTYPE = c()
for(i in 1:N)
{
  pheno_name = PHENOTYPES[i]
  if (i%%10==0)
  {
    print(paste0(">>> Deal with phenotype ",pheno_name," (",i,"/",N,"). Nb collected = ",length(COLLECTED_PHENOTYPE)))
  }
  file                = paste0("/scratch/project_XXXXXXX/Tribolium_castaneum_eQTL/rds/",DATASET,"_",pheno_name,".rds")
  d                   = readRDS(file)
  d                   = d[!is.na(d$p_wald),]
  d$gc                = calculate_gc(d$p_wald)
  d$fdr               = calculate_fdr(d$gc)
  THRESHOLD           = 0.05
  d_select            = d[d$fdr<THRESHOLD,]
  COLLECTED_EQTL      = rbind(COLLECTED_EQTL, d_select)
  COLLECTED_PHENOTYPE = c(COLLECTED_PHENOTYPE, rep(pheno_name, dim(d_select)[1]))
}
COLLECTED_EQTL$phenotype = COLLECTED_PHENOTYPE

#--------------------------------#
# 4) Save resulting datasets     #
#--------------------------------#
fileRDS = paste0("/scratch/project_XXXXXXX/Tribolium_castaneum_eQTL/",DATASET,"_significant.rds")
upload_to_allas(fileRDS, "ecoevodyn_tribolium_tcas3_30")
saveRDS(COLLECTED_EQTL, file=fileRDS)

fileCSV = paste0("/scratch/project_XXXXXXX/Tribolium_castaneum_eQTL/",DATASET,"_significant.csv")
upload_to_allas(fileCSV, "ecoevodyn_tribolium_tcas3_30")
write.table(COLLECTED_EQTL, file=fileCSV, row.names=F, col.names=T, quote=F, sep="\t")

