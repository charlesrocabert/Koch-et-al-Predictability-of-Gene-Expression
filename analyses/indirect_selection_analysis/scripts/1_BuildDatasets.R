#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 1_BuildDatasets.R
# -----------------
# Build SNP, gene and eQTL datasets prior to analyses.
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")
library('rstudioapi')

### Mean wrapper function ###
my_mean <- function( x )
{
  return(mean(x, na.rm=T))
}

### Sum wrapper function ###
my_sum <- function( x )
{
  return(sum(x, na.rm=T))
}

### Load eQTLs ###
load_eQTLs <- function( phenotype )
{
  #--------------------------#
  # 1) Load eQTLs data       #
  #--------------------------#
  EQTLs    = readRDS(paste0("./data/tribolium_eqtl/significant/HD_G1_Tcas3.30_imputed_",phenotype,"_significant.rds"))
  EQTLs$ID = EQTLs$rs
  #--------------------------#
  # 2) Load SNP annotation   #
  #--------------------------#
  ANNOTATION        = read.table("./data/tribolium_snp/snp_table_ALL_Tcas3.30_raw_SNP.csv", h=T, sep="\t")
  ANNOTATION$gene   = ANNOTATION$Feature_id
  ANNOTATION_SELECT = select(ANNOTATION, c("ID", "Annotation", "Putative_impact", "gene"))
  #--------------------------#
  # 3) Select and merge data #
  #--------------------------#
  eQTLs = merge(EQTLs, ANNOTATION_SELECT, by="ID")
  #--------------------------#
  # 4) Return data           #
  #--------------------------#
  rm(ANNOTATION)
  return(eQTLs)
}

### Build the SNP dataset ###
build_SNP_dataset <- function( eQTLs )
{
  #-------------------#
  # 1) Load datasets  #
  #-------------------#
  
  ### Load allelic frequency changes (AFCs) data ###
  AFC_HD = readRDS(paste0("./data/tribolium_afc/AFC_HD.rds"))
  AFC_HD = AFC_HD[!duplicated(AFC_HD$ID),]
  AFC_HD = AFC_HD[,c("ID","POS","CHROM")]
  
  ### Load SNP parallelism information ###
  L1_SIGNIF  = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_L1.rds")
  L2_SIGNIF  = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_L2.rds")
  L3_SIGNIF  = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_L3.rds")
  L5_SIGNIF  = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_L5.rds")
  L6_SIGNIF  = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_L6.rds")
  Mx1_SIGNIF = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_Mx1.rds")
  Mx2_SIGNIF = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_Mx2.rds")
  ########################################
  L1_SIGNIF  = as.data.frame(L1_SIGNIF$alpha05)
  L2_SIGNIF  = as.data.frame(L2_SIGNIF$alpha05)
  L3_SIGNIF  = as.data.frame(L3_SIGNIF$alpha05)
  L5_SIGNIF  = as.data.frame(L5_SIGNIF$alpha05)
  L6_SIGNIF  = as.data.frame(L6_SIGNIF$alpha05)
  Mx1_SIGNIF = as.data.frame(Mx1_SIGNIF$alpha05)
  Mx2_SIGNIF = as.data.frame(Mx2_SIGNIF$alpha05)
  ########################################
  names(L1_SIGNIF)  = c("index", "pvalue")
  names(L2_SIGNIF)  = c("index", "pvalue")
  names(L3_SIGNIF)  = c("index", "pvalue")
  names(L5_SIGNIF)  = c("index", "pvalue")
  names(L6_SIGNIF)  = c("index", "pvalue")
  names(Mx1_SIGNIF) = c("index", "pvalue")
  names(Mx2_SIGNIF) = c("index", "pvalue")
  ########################################
  L1_SIGNIF$index   = as.numeric(L1_SIGNIF$index)
  L1_SIGNIF$pvalue  = as.numeric(L1_SIGNIF$pvalue)
  L2_SIGNIF$index   = as.numeric(L2_SIGNIF$index)
  L2_SIGNIF$pvalue  = as.numeric(L2_SIGNIF$pvalue)
  L3_SIGNIF$index   = as.numeric(L3_SIGNIF$index)
  L3_SIGNIF$pvalue  = as.numeric(L3_SIGNIF$pvalue)
  L5_SIGNIF$index   = as.numeric(L5_SIGNIF$index)
  L5_SIGNIF$pvalue  = as.numeric(L5_SIGNIF$pvalue)
  L6_SIGNIF$index   = as.numeric(L6_SIGNIF$index)
  L6_SIGNIF$pvalue  = as.numeric(L6_SIGNIF$pvalue)
  Mx1_SIGNIF$index  = as.numeric(Mx1_SIGNIF$index)
  Mx1_SIGNIF$pvalue = as.numeric(Mx1_SIGNIF$pvalue)
  Mx2_SIGNIF$index  = as.numeric(Mx2_SIGNIF$index)
  Mx2_SIGNIF$pvalue = as.numeric(Mx2_SIGNIF$pvalue)
  
  ### Load SNP annotation ###
  ANNOTATION        = read.table("./data/tribolium_snp/snp_table_ALL_Tcas3.30_raw_SNP.csv", h=T, sep="\t")
  ANNOTATION$gene   = ANNOTATION$Feature_id
  ANNOTATION_SELECT = select(ANNOTATION, c("ID", "Annotation", "Putative_impact", "gene"))
  
  ### Load list of genes under selection ###
  
  # 1) Genes with a significant s_g AND Genes with the highest abs(s_g)
  # -------------------------------------------------------------------
  GENE_FITNESS_COR_DATASET          = read.table("./data/experiment_data/GeneticCor_Fitness_HD.csv", sep=";", h=T, dec=",")
  GENE_FITNESS_COR_DATASET$Mean_abs = abs(GENE_FITNESS_COR_DATASET$Mean)
  GENE_FITNESS_COR_DATASET          = GENE_FITNESS_COR_DATASET[order(GENE_FITNESS_COR_DATASET$Mean_abs, decreasing=T),]
  significant_sg_pool               = filter(GENE_FITNESS_COR_DATASET, sig>0)$gene
  highest_sg_pool                   = GENE_FITNESS_COR_DATASET$gene[1:106]
  
  # 2) Genes showing an evolutionary response (Koch & Guillaume, 2020)
  # ------------------------------------------------------------------
  EVOL_RESPONSE_DATASET           = read.table("./data/experiment_data/Evol-Diff-Expr_HD_G22.txt", sep="\t", h=T)
  EVOL_RESPONSE_DATASET$abs_logFC = abs(EVOL_RESPONSE_DATASET$logFC)
  EVOL_RESPONSE_DATASET           = EVOL_RESPONSE_DATASET[order(EVOL_RESPONSE_DATASET$abs_logFC, decreasing=T),]
  significant_de_pool             = rownames(filter(EVOL_RESPONSE_DATASET, DE!=0))
  highest_de_pool                 = rownames(EVOL_RESPONSE_DATASET[1:135,])
  
  # 3) Create a common group for all previous sets
  # ----------------------------------------------
  common_pool = c(significant_sg_pool, highest_sg_pool, significant_de_pool, highest_de_pool)
  
  # 4) Genes in expression modules correlated to fitness
  # ----------------------------------------------------
  EXPRESSION_MODULES      = read.table("./data/tribolium_modules/HD_G1_Tcas3.30_EXPRESSION_significant_modules.txt", sep="\t", h=T)
  expression_modules_pool = EXPRESSION_MODULES$gene
  
  # 5) WGCNA hub genes
  # ------------------
  gene_connectivity = read.table("./data/tribolium_modules/results_eva/Connectivity_HD.txt", sep="\t", h=T)
  gene_connectivity = gene_connectivity[order(gene_connectivity$Module, gene_connectivity$kWithin, decreasing=T),]
  hub_genes         = c()
  module_genes      = c()
  for(module in unique(gene_connectivity$Module))
  {
    N            = length(rownames(filter(gene_connectivity, Module==module)))
    hub_genes    = c(hub_genes, rownames(filter(gene_connectivity, Module==module))[1:10])
    module_genes = c(module_genes, rownames(filter(gene_connectivity, Module==module))[11:N])
  }
  
  ### Load phenotype positions ###
  PHENOTYPE_POSITIONS        = read.table("./data/tribolium_eqtl/gene_pos_Tcas3.30.csv", sep="\t", h=T)
  names(PHENOTYPE_POSITIONS) = c("gene", "gene_chr", "gene_start", "gene_end")
  PHENOTYPE_POSITIONS        = PHENOTYPE_POSITIONS[,c("gene", "gene_start", "gene_end")]
  
  ### Load fitness covariance results ###
  FITNESS_DATA           = read.table("./data/experiment_data/Selection-Cov-Beta-HD.csv", h=T, sep=",")
  rownames(FITNESS_DATA) = FITNESS_DATA$gene
  
  #-------------------#
  # 2) Merge datasets #
  #-------------------#
  
  ### Create the final dataset ###
  DATA = data.frame(AFC_HD)
  
  ### Merge parallelism information ###
  DATA$L1_shift            = as.numeric(DATA$ID%in%DATA$ID[L1_SIGNIF$index])
  DATA$L2_shift            = as.numeric(DATA$ID%in%DATA$ID[L2_SIGNIF$index])
  DATA$L3_shift            = as.numeric(DATA$ID%in%DATA$ID[L3_SIGNIF$index])
  DATA$L5_shift            = as.numeric(DATA$ID%in%DATA$ID[L5_SIGNIF$index])
  DATA$L6_shift            = as.numeric(DATA$ID%in%DATA$ID[L6_SIGNIF$index])
  DATA$Mx1_shift           = as.numeric(DATA$ID%in%DATA$ID[Mx1_SIGNIF$index])
  DATA$Mx2_shift           = as.numeric(DATA$ID%in%DATA$ID[Mx2_SIGNIF$index])
  DATA$Parallelism         = DATA$L1_shift+DATA$L2_shift+DATA$L3_shift+DATA$L5_shift+DATA$L6_shift+DATA$Mx1_shift+DATA$Mx2_shift
  DATA$isShifting          = as.numeric(DATA$Parallelism>0)
  DATA$isPartiallyParallel = as.numeric(DATA$Parallelism>1 & DATA$Parallelism<5)
  DATA$isFullyParallel     = as.numeric(DATA$Parallelism>4)
  
  ### Build parallel categories ###
  DATA$Parallel_category = as.character(DATA$Parallelism)
  DATA$Parallel_category[DATA$Parallel_category=="0"] = "No shift"
  DATA$Parallel_category[DATA$Parallel_category=="1"] = "Single shift"
  DATA$Parallel_category[DATA$Parallel_category=="2"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="3"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="4"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="5"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="6"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="7"] = "Partial/full parallelism"
  DATA$Parallel_category = factor(DATA$Parallel_category, levels=c("No shift","Single shift","Partial/full parallelism"))
  
  ### Merge annotation ###
  DATA = merge(DATA, ANNOTATION_SELECT, by="ID")
  
  ### Merge with evolutionary response ###
  DATA$logFC     = EVOL_RESPONSE_DATASET[DATA$gene, "logFC"]
  DATA$abs_logFC = EVOL_RESPONSE_DATASET[DATA$gene, "abs_logFC"]
  
  ### Merge with genes under selection ###
  DATA$significant_sg_pool = as.numeric(DATA$gene%in%significant_sg_pool)
  DATA$highest_sg_pool     = as.numeric(DATA$gene%in%highest_sg_pool)
  DATA$significant_de_pool = as.numeric(DATA$gene%in%significant_de_pool)
  DATA$highest_de_pool     = as.numeric(DATA$gene%in%highest_de_pool)
  DATA$common_pool         = as.numeric(DATA$gene%in%common_pool)
  DATA$hub_gene            = as.numeric(DATA$gene%in%hub_genes)
  DATA$module_gene         = as.numeric(DATA$gene%in%module_genes)
  
  ### Merge with gene positions ###
  DATA = merge(DATA, PHENOTYPE_POSITIONS, by="gene")
  
  ### Merge with fitness data ###
  # COV_g  = direct + undirect (net)
  # beta_g = direct only
  #
  DATA$total_selection      = FITNESS_DATA[DATA$gene,"Cov.Mean"]
  DATA$direct_selection     = FITNESS_DATA[DATA$gene,"bgHD.Mean"]
  DATA$diff_selection       = DATA$total_selection-DATA$direct_selection
  DATA$log_total_selection  = log10(abs(DATA$total_selection))
  DATA$log_direct_selection = log10(abs(DATA$direct_selection))
  
  ### Which SNPs are eQTLs ###
  DATA$is_eQTL = as.numeric(DATA$ID %in% eQTLs$ID)
  
  ### Clear memory ###
  rm(AFC_HD)
  rm(L1_SIGNIF)
  rm(L2_SIGNIF)
  rm(L3_SIGNIF)
  rm(L5_SIGNIF)
  rm(L6_SIGNIF)
  rm(Mx1_SIGNIF)
  rm(Mx2_SIGNIF)
  rm(ANNOTATION)
  rm(ANNOTATION_SELECT)
  rm(EVOL_RESPONSE_DATASET)
  rm(GENE_FITNESS_COR_DATASET)
  rm(PHENOTYPE_POSITIONS)
  rm(FITNESS_DATA)
  
  ### Return the dataset ###
  rownames(DATA) = DATA$ID
  return(DATA)
}

### Build the gene dataset ###
build_gene_dataset <- function( SNP_dataset, eQTLs )
{
  to_remove = c("ID", "POS", "L1_shift", "L2_shift", "L3_shift", "L5_shift", "L6_shift", "Mx1_shift", "Mx2_shift", "is_eQTL")
  ### Collapse data to genes ###
  DATA                     = SNP_dataset[!duplicated(SNP_dataset$gene),]
  DATA                     = DATA[,-which(names(DATA)%in%to_remove)]
  DATA$Parallelism         = tapply(SNP_dataset$Parallelism, SNP_dataset$gene, max)[DATA$gene]
  DATA$isShifting          = tapply(SNP_dataset$isShifting, SNP_dataset$gene, sum)[DATA$gene]
  DATA$isShifting          = as.numeric(DATA$isShifting>0)
  DATA$isPartiallyParallel = tapply(SNP_dataset$isPartiallyParallel, SNP_dataset$gene, sum)[DATA$gene]
  DATA$isPartiallyParallel = as.numeric(DATA$isPartiallyParallel>0)
  DATA$isFullyParallel     = tapply(SNP_dataset$isFullyParallel, SNP_dataset$gene, sum)[DATA$gene]
  DATA$isFullyParallel     = as.numeric(DATA$isFullyParallel>0)
  
  ### Build parallel categories ###
  DATA$Parallel_category = as.character(DATA$Parallelism)
  DATA$Parallel_category[DATA$Parallel_category=="0"] = "No shift"
  DATA$Parallel_category[DATA$Parallel_category=="1"] = "Single shift"
  DATA$Parallel_category[DATA$Parallel_category=="2"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="3"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="4"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="5"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="6"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="7"] = "Partial/full parallelism"
  DATA$Parallel_category = factor(DATA$Parallel_category, levels=c("No shift","Single shift","Partial/full parallelism"))
  
  ### Which genes are hub genes ###
  DATA$hub_gene_category = as.character(DATA$hub_gene)
  DATA$hub_gene_category[DATA$hub_gene_category=="0"] = "Other"
  DATA$hub_gene_category[DATA$hub_gene_category=="1"] = "Hub gene"
  DATA$hub_gene_category = factor(DATA$hub_gene_category, levels=c("Other","Hub gene"))
  
  ### Return the dataset ###
  rownames(DATA) = DATA$gene
  return(DATA)
}

### Build the eQTL dataset ###
build_eQTL_dataset <- function( SNP_dataset, gene_dataset, eQTLs )
{
  ### Build the initial eQTL dataset ###
  DATA                           = data.frame(eQTLs)
  to_keep                        = names(DATA)[!names(DATA)%in%c("Annotation","Putative_impact","gene")]
  DATA                           = DATA[,to_keep]
  DATA                           = merge(DATA, SNP_dataset, by="ID")
  DATA$total_selection           = gene_dataset[DATA$phenotype, "total_selection"]
  DATA$direct_selection          = gene_dataset[DATA$phenotype, "direct_selection"]
  DATA$weighted_total_selection  = DATA$beta*((DATA$total_selection))
  DATA$weighted_direct_selection = DATA$beta*((DATA$direct_selection))
  
  ### Calculate additional statistics ###
  mean_total_selection           = tapply(DATA$total_selection, DATA$ID, my_mean)
  mean_direct_selection          = tapply(DATA$direct_selection, DATA$ID, my_mean)
  weighted_total_selection       = tapply(DATA$weighted_total_selection, DATA$ID, my_sum)
  weighted_direct_selection      = tapply(DATA$weighted_direct_selection, DATA$ID, my_sum)
  DATA$Pleiotropy                = as.vector(table(DATA$ID)[DATA$ID])
  DATA$mean_total_selection      = mean_total_selection[DATA$ID]
  DATA$mean_direct_selection     = mean_direct_selection[DATA$ID]
  DATA$weighted_total_selection  = weighted_total_selection[DATA$ID]
  DATA$weighted_direct_selection = weighted_direct_selection[DATA$ID]
  DATA$log_total_selection       = log10(abs(DATA$weighted_total_selection))
  DATA$log_direct_selection      = log10(abs(DATA$weighted_direct_selection))
  
  ### Build pleiotropy classes ###
  DATA$Pleiotropy_category = DATA$Pleiotropy
  DATA$Pleiotropy_category[DATA$Pleiotropy_category>1 & DATA$Pleiotropy_category<5] = 2
  DATA$Pleiotropy_category[DATA$Pleiotropy_category>=5] = 3
  DATA$Pleiotropy_category = as.character(DATA$Pleiotropy_category)
  DATA$Pleiotropy_category[DATA$Pleiotropy_category=="1"] = "No pleiotropy"
  DATA$Pleiotropy_category[DATA$Pleiotropy_category=="2"] = "Low pleiotropy"
  DATA$Pleiotropy_category[DATA$Pleiotropy_category=="3"] = "High pleiotropy"
  DATA$Pleiotropy_category = factor(DATA$Pleiotropy_category, levels=c("No pleiotropy","Low pleiotropy","High pleiotropy"))
  
  ### Build parallel classes ###
  DATA$Parallel_category = as.character(DATA$Parallelism)
  DATA$Parallel_category[DATA$Parallel_category=="0"] = "No shift"
  DATA$Parallel_category[DATA$Parallel_category=="1"] = "Single shift"
  DATA$Parallel_category[DATA$Parallel_category=="2"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="3"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="4"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="5"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="6"] = "Partial/full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="7"] = "Partial/full parallelism"
  DATA$Parallel_category = factor(DATA$Parallel_category, levels=c("No shift","Single shift","Partial/full parallelism"))
  
  ### Clear memory ###
  rm(mean_total_selection)
  rm(mean_direct_selection)
  rm(weighted_total_selection)
  rm(weighted_direct_selection)
  
  ### Return the dataset ###
  #rownames(DATA) = DATA$ID
  return(DATA)
}

### Build the eQTL carrier dataset ###
build_eQTL_carrier_dataset <- function( eQTL_dataset, gene_dataset )
{
  ### Calculate pleiotropy at the gene level ###
  Pleiotropy = c()
  for (g in unique(eQTL_dataset$gene))
  {
    val        = length(unique(eQTL_dataset[eQTL_dataset$gene==g, "phenotype"]))
    Pleiotropy = c(Pleiotropy, val)
  }
  names(Pleiotropy) = unique(eQTL_dataset$gene)
  
  ### Merge datasets ###
  DATA           = data.frame(names(Pleiotropy), Pleiotropy)
  names(DATA)    = c("gene", "Pleiotropy")
  rownames(DATA) = DATA$gene
  
  ### Return the dataset ###
  return(DATA)
}

### Build the eQTL phenotype dataset ###
build_eQTL_phenotype_dataset <- function( SNP_dataset, gene_dataset, eQTLs )
{
  ### Calculate phenotype connectivity at the gene level ###
  Connectivity = c()
  for (p in unique(eQTLs$phenotype))
  {
    val          = length(unique(eQTLs[eQTLs$phenotype==p, "gene"]))
    Connectivity = c(Connectivity, val)
  }
  names(Connectivity) = unique(eQTLs$phenotype)
  
  ### Merge datasets ###
  DATA           = data.frame(names(Connectivity), Connectivity)
  names(DATA)    = c("gene", "Connectivity")
  rownames(DATA) = DATA$gene
  
  ### Return the dataset ###
  return(DATA)
}

### Add eQTL pleiotropy and eQTL phenotype connectivity categories ###
add_pleiotropy_connectivity <- function ( gene_dataset, eQTL_carrier_dataset, eQTL_phenotype_dataset )
{
  ### What is the pleiotropy of each gene? (0 if the gene is not a hub gene) ###
  gene_dataset$Pleiotropy = eQTL_carrier_dataset[gene_dataset$gene, "Pleiotropy"]
  gene_dataset$Pleiotropy[is.na(gene_dataset$Pleiotropy)] = 0
  gene_dataset$eQTL_carrier = as.numeric(gene_dataset$Pleiotropy > 0)
  
  gene_dataset$Pleiotropy_category = gene_dataset$Pleiotropy
  gene_dataset$Pleiotropy_category[gene_dataset$Pleiotropy_category>1] = 2
  gene_dataset$Pleiotropy_category = as.character(gene_dataset$Pleiotropy_category)
  gene_dataset$Pleiotropy_category[gene_dataset$Pleiotropy_category=="0"] = "Non-eQTL carrier"
  gene_dataset$Pleiotropy_category[gene_dataset$Pleiotropy_category=="1"] = "Single-phenotype eQTL carrier"
  gene_dataset$Pleiotropy_category[gene_dataset$Pleiotropy_category=="2"] = "Pleiotropic eQTL carrier"
  gene_dataset$Pleiotropy_category = factor(gene_dataset$Pleiotropy_category, levels=c("Non-eQTL carrier", "Single-phenotype eQTL carrier","Pleiotropic eQTL carrier"))
  
  ### What is the connectivity of each gene? (0 if the gene is not a eQTL phenotype) ###
  gene_dataset$Connectivity = eQTL_phenotype_dataset[gene_dataset$gene, "Connectivity"]
  gene_dataset$Connectivity[is.na(gene_dataset$Connectivity)] = 0
  gene_dataset$eQTL_phenotype = as.numeric(gene_dataset$Connectivity > 0)
  
  ### Build connectivity classes ###
  gene_dataset$Connectivity_category = gene_dataset$Connectivity
  gene_dataset$Connectivity_category[gene_dataset$Connectivity_category>0 & gene_dataset$Connectivity_category<5] = 1
  gene_dataset$Connectivity_category[gene_dataset$Connectivity_category>4] = 2
  gene_dataset$Connectivity_category = as.character(gene_dataset$Connectivity_category)
  gene_dataset$Connectivity_category[gene_dataset$Connectivity_category=="0"] = "Non-eQTL phenotype"
  gene_dataset$Connectivity_category[gene_dataset$Connectivity_category=="1"] = "Low connectivity"
  gene_dataset$Connectivity_category[gene_dataset$Connectivity_category=="2"] = "High connectivity"
  gene_dataset$Connectivity_category = factor(gene_dataset$Connectivity_category, levels=c("Non-eQTL phenotype","Low connectivity","High connectivity"))
  
  ### Return the updated gene dataset ###
  return(gene_dataset)
}


##################
#      MAIN      #
##################

directory = dirname(getActiveDocumentContext()$path)
setwd(directory)
setwd("../../..")

#----------------------------------------------#
# 1) Build the new datasets                    #
#----------------------------------------------#
eQTLs                  = load_eQTLs("EXPRESSION")
SNP_dataset            = build_SNP_dataset(eQTLs)
gene_dataset           = build_gene_dataset(SNP_dataset, eQTLs)
eQTL_dataset           = build_eQTL_dataset(SNP_dataset, gene_dataset, eQTLs)
eQTL_carrier_dataset   = build_eQTL_carrier_dataset(eQTL_dataset, gene_dataset)
eQTL_phenotype_dataset = build_eQTL_phenotype_dataset(SNP_dataset, gene_dataset, eQTLs)

#----------------------------------------------#
# 2) Complete the information                  #
#----------------------------------------------#
gene_dataset = add_pleiotropy_connectivity(gene_dataset, eQTL_carrier_dataset, eQTL_phenotype_dataset)

#----------------------------------------------#
# 3) Save the datasets                         #
#----------------------------------------------#
saveRDS(SNP_dataset, "./analyses/indirect_selection_analysis/data/SNP_dataset.rds")
saveRDS(gene_dataset, "./analyses/indirect_selection_analysis/data/gene_dataset.rds")
saveRDS(eQTL_phenotype_dataset, "./analyses/indirect_selection_analysis/data/eQTL_phenotype_dataset.rds")
saveRDS(eQTL_carrier_dataset, "./analyses/indirect_selection_analysis/data/eQTL_carrier_dataset.rds")

#----------------------------------------------#
# 4) Save various gene lists for GO enrichment #
#----------------------------------------------#
significant_sg_genes = unique(gene_dataset[gene_dataset$significant_sg_pool==1, "gene"])
highest_sg_genes     = unique(gene_dataset[gene_dataset$highest_sg_pool==1, "gene"])
significant_de_genes = unique(gene_dataset[gene_dataset$significant_de_pool==1, "gene"])
highest_de_genes     = unique(gene_dataset[gene_dataset$highest_de_pool==1, "gene"])
common_pool_genes    = unique(gene_dataset[gene_dataset$common_pool==1, "gene"])
parallel_genes       = unique(gene_dataset[gene_dataset$isFullyParallel==1, "gene"])
eQTL_carriers        = unique(gene_dataset[gene_dataset$eQTL_carrier==1, "gene"])
eQTL_phenotypes      = unique(gene_dataset[gene_dataset$eQTL_phenotype==1, "gene"])
hub_genes            = unique(gene_dataset[gene_dataset$hub_gene==1, "gene"])

write.table(significant_sg_genes, "./analyses/indirect_selection_analysis/data/gene_sets/significant_sg_genes_list.txt", quote=F, row.names=F, col.names=F)
write.table(highest_sg_genes, "./analyses/indirect_selection_analysis/data/gene_sets/highest_sg_genes_list.txt", quote=F, row.names=F, col.names=F)
write.table(significant_de_genes, "./analyses/indirect_selection_analysis/data/gene_sets/significant_de_genes_list.txt", quote=F, row.names=F, col.names=F)
write.table(highest_de_genes, "./analyses/indirect_selection_analysis/data/gene_sets/highest_de_genes_list.txt", quote=F, row.names=F, col.names=F)
write.table(common_pool_genes, "./analyses/indirect_selection_analysis/data/gene_sets/common_pool_genes_list.txt", quote=F, row.names=F, col.names=F)
write.table(parallel_genes, "./analyses/indirect_selection_analysis/data/gene_sets/parallel_genes_list.txt", quote=F, row.names=F, col.names=F)
write.table(eQTL_carriers, "./analyses/indirect_selection_analysis/data/gene_sets/eQTL_carriers_list.txt", quote=F, row.names=F, col.names=F)
write.table(eQTL_phenotypes, "./analyses/indirect_selection_analysis/data/gene_sets/eQTL_phenotypes_list.txt", quote=F, row.names=F, col.names=F)
write.table(hub_genes, "./analyses/indirect_selection_analysis/data/gene_sets/hub_genes_list.txt", quote=F, row.names=F, col.names=F)

