#!/usr/bin/env Rscript

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# 2_ShowStatistics.R
# ------------------
# Display some statistics for the manuscript.
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

library("tidyverse")
library("rstudioapi")
library("ggcorrplot")

### Return the stars corresponding to pvalues ###
get_pvalue_stars <- function( pvalue )
{
  star = "ns"
  if (pvalue < 0.0001)
  {
    star = "****"
  }
  if (pvalue >= 0.0001 & pvalue < 0.001)
  {
    star = "***"
  }
  if (pvalue >= 0.001 & pvalue < 0.01)
  {
    star = "**"
  }
  if (pvalue >= 0.01 & pvalue < 0.05)
  {
    star = "*"
  }
  if (pvalue >= 0.05 & pvalue < 0.1)
  {
    star = "."
  }
  if (pvalue >= 0.1)
  {
    star = "ns"
  }
  return(star)
}

### Aggregator function ###
aggregator_function <- function( x, aggregator )
{
  if (aggregator == "sum")
  {
    return(sum(x, na.rm=T))
  }
  if (aggregator == "mean")
  {
    return(mean(x, na.rm=T))
  }
  if (aggregator == "median")
  {
    return(median(x, na.rm=T))
  }
}

### Run one permutations test ###
run_permutations_test <- function( dataset, target_variable, source_variable, source_list, nb_reps, aggregator )
{
  x              = dataset[,target_variable]
  pos            = which(dataset[,source_variable]%in%source_list)
  original_score = aggregator_function(x[pos], aggregator)
  permutations   = t(sapply(1:nb_reps,function(i){
    score = aggregator_function(sample(x)[pos], aggregator)
    return(score)
  }))
  qval = sum(permutations<original_score)/nb_reps
  return(1-qval)
}


##################
#      MAIN      #
##################

directory = dirname(getActiveDocumentContext()$path)
setwd(directory)
setwd("../../..")

#-------------------------------------------------------#
# 1) Load the data                                      #
#-------------------------------------------------------#
SNP_dataset            = readRDS("./analyses/indirect_selection_analysis/data/SNP_dataset.rds")
gene_dataset           = readRDS("./analyses/indirect_selection_analysis/data/gene_dataset.rds")
eQTL_dataset           = readRDS("./analyses/indirect_selection_analysis/data/eQTL_dataset.rds")
eQTL_carrier_dataset   = readRDS("./analyses/indirect_selection_analysis/data/eQTL_carrier_dataset.rds")
eQTL_phenotype_dataset = readRDS("./analyses/indirect_selection_analysis/data/eQTL_phenotype_dataset.rds")

#-------------------------------------------------------#
# 2) Display headers for supp. mat. annotation          #
#-------------------------------------------------------#
line = ""
for (n in names(gene_dataset))
{
  line = paste0(line, "- ", n, "\n")
}
cat(line)

#-------------------------------------------------------#
# 3) Display some general information                   #
#-------------------------------------------------------#
cat(paste0("• Total number of SNPs: ", length(unique(SNP_dataset$ID)),"\n",
          "• Total number of genes: ", length(unique(gene_dataset$gene)),"\n",
          "• Total number of eQTLs: ", length(unique(eQTL_dataset$ID)),"\n",
          "• Total number of phenotypes: ", length(unique(eQTL_dataset$phenotype)),"\n",
          "• Total number of eQTLs carrier genes: ", length(unique(eQTL_carrier_dataset$gene)),"\n"))
####################
cat(paste0("• SNPs shifting in at least one line: ", length(unique(filter(SNP_dataset, Parallelism>0)$ID)),"\n",
           "• SNPs with parallelism = 1: ", length(unique(filter(SNP_dataset, Parallelism==1)$ID)),"\n",
           "• SNPs with partial parallelism: ", length(unique(filter(SNP_dataset, Parallelism>1 & Parallelism<5)$ID)),"\n",
           "• SNPs with high parallelism (>5): ", length(unique(filter(SNP_dataset, Parallelism>4)$ID)),"\n"))
####################
cat(paste0("• Genes shifting in at least one line: ", length(unique(filter(gene_dataset, Parallelism>0)$ID)),"\n",
           "• Genes with parallelism = 1: ", length(unique(filter(gene_dataset, Parallelism==1)$ID)),"\n",
           "• Genes with partial parallelism: ", length(unique(filter(gene_dataset, Parallelism>1 & Parallelism<5)$ID)),"\n",
           "• Genes with high parallelism (>5): ", length(unique(filter(gene_dataset, Parallelism>4)$ID)),"\n"))

#-------------------------------------------------------#
# 4) Are hub genes enriched in parallel genes?          #
#-------------------------------------------------------#

### Permutation test ###
source_list = filter(gene_dataset, hub_gene==1)$gene
run_permutations_test(gene_dataset, "isFullyParallel", "gene", source_list, 10000, "sum")

### Chi-square test ###
gene_dataset$ParallelismDegree = gene_dataset$Parallelism
gene_dataset$ParallelismDegree[gene_dataset$ParallelismDegree>1] = 2
X = table(gene_dataset$ParallelismDegree, gene_dataset$hub_gene)
chisq.test(X)

#-------------------------------------------------------#
# 5) Are hub genes enriched in eQTLs?                   #
#-------------------------------------------------------#
X = table(SNP_dataset$is_eQTL, SNP_dataset$hub_gene)
chisq.test(X)

#-------------------------------------------------------#
# 6) Are significant sg genes enriched in parallel AFC? #
#-------------------------------------------------------#
source_list = filter(SNP_dataset, significant_sg_pool==1)$ID
run_permutations_test(SNP_dataset, "isFullyParallel", "ID", source_list, 10000, "sum")

#-------------------------------------------------------#
# 7) Are eQTL carriers more parallel?                   #
#-------------------------------------------------------#
X = table(gene_dataset$isFullyParallel, gene_dataset$eQTL_carrier)
chisq.test(X)

