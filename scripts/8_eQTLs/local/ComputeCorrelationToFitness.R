#!/usr/bin/env Rscript

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# ComputeCorrelationToFitness.R
# -----------------------------
# Compute correlation of each phenotype to fitness.
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

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

#--------------------------------#
# 1) Read command line arguments #
#--------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
POPULATION      = args[2]
VERSION         = args[3]
PHENOTYPE       = args[4]
setwd(REPOSITORY_PATH)

#--------------------------------#
# 2) Load gene expression data   #
#--------------------------------#
expr           = read.table(paste0("./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_",PHENOTYPE,".txt"), h=T, sep="\t", check.names=F)
gene_id        = expr[,1]
expr           = expr[,2:dim(expr)[2]]
samples        = colnames(expr)
rownames(expr) = gene_id

#--------------------------------#
# 3) Load fitness data           #
#--------------------------------#
fitness = read.table(paste0("./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_FITNESS.txt"), h=T, sep="\t", check.names=F)
fitness = fitness[samples,2]

#--------------------------------#
# 4) Compute correlations        #
#--------------------------------#
ESTIMATES = c()
PVALUES   = c()
N         = length(gene_id)
for (i in 1:N)
{
  gene = gene_id[i]
  if (i%%10==0)
  {
    print(paste0(">>> Deal with phenotype ",gene," (",i,"/",N,")"))
  }
  estimate  = cor.test(t(expr[gene,]), fitness)$estimate
  pvalue    = cor.test(t(expr[gene,]), fitness)$p.value
  ESTIMATES = c(ESTIMATES, estimate)
  PVALUES   = c(PVALUES, pvalue)
}
dataset = data.frame(gene_id, ESTIMATES, PVALUES)
names(dataset) = c("phenotype", "estimate", "pvalue")

#--------------------------------#
# 5) Save the dataset            #
#--------------------------------#
filename = paste0(POPULATION,"_",VERSION,"_",PHENOTYPE,"_fitness_cor.rds")
saveRDS(dataset, file=paste0("./data/tribolium_eqtl/fitness_cor/",filename))

