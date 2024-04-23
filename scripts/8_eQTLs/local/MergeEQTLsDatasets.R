#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# MergeEQTLsDatasets.R
# --------------------
# Merge eQTLs data with other information (SNP annotation and fitness).
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

#--------------------------------------------#
# 1) Read command line arguments             #
#--------------------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<6)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
POPULATION      = args[2]
VERSION         = args[3]
IMPUTATION      = args[4]
if (IMPUTATION=="-")
{
  IMPUTATION = ""
}
PHENOTYPE   = args[5]
ENVIRONMENT = args[6]
setwd(REPOSITORY_PATH)

#--------------------------------------------#
# 2) Load significant eQTLs                  #
#--------------------------------------------#
EQTL            = readRDS(paste0("./data/tribolium_eqtl/significant/",POPULATION,"_",VERSION,IMPUTATION,"_",PHENOTYPE,"_significant.rds"))
names(EQTL)[2]  = "ID"
names(EQTL)[3]  = "pos"
names(EQTL)[12] = "p"
MERGED          = data.frame(EQTL)

rm(EQTL)

#--------------------------------------------#
# 3) Load phenotype fitness correlation data #
#--------------------------------------------#
FITNESS           = readRDS(paste0("./data/tribolium_eqtl/fitness_cor/",POPULATION,"_",VERSION,"_",PHENOTYPE,"_fitness_cor.rds"))
rownames(FITNESS) = FITNESS$phenotype
names(FITNESS)    = c("phenotype", "fitness_cor_estimate", "fitness_cor_pvalue")
MERGED            = cbind(MERGED, FITNESS[MERGED$phenotype,c("fitness_cor_estimate","fitness_cor_pvalue")])

rm(FITNESS)

#--------------------------------------------#
# 4) Extract G1 SNPs annotation data         #
#--------------------------------------------#

### Load SNP data ###
SNP_G1 = read.table(paste0("./data/tribolium_snp/snp_table_",ENVIRONMENT,"_G1_",VERSION,"_eQTL",IMPUTATION,".csv"), h=T, sep="\t", check.names=F)

### Extract pertinent information ###
SNP_G1 = SNP_G1[,c("ID", "Annotation", "Putative_impact", "Gene_name", "Gene_id", "Feature_type", "Feature_id", "Transcript_biotype", "Rank_total", "HGVS_c", "HGVS_p", "cDNA_position", "CDS_position", "Protein_position", "Distance_to_feature", "Errors")]

### Rename certain columns ###
names(SNP_G1) = c("ID", "Annotation", "Putative_impact", "Gene_name", "Gene_id", "Feature_type", "Feature_id", "Transcript_biotype", "Rank_total", "HGVS_c", "HGVS_p", "cDNA_position", "CDS_position", "Protein_position", "Distance_to_feature", "Errors")

### Merge with main dataset ###
MERGED = merge(MERGED, SNP_G1, by="ID")

rm(SNP_G1)

#--------------------------------------------#
# 5) Determine cis-trans status              #
#--------------------------------------------#

### Load gene positions ###
gene_pos = read.table(paste0("./data/tribolium_eqtl/gene_pos_",VERSION,".csv"), h=T, sep="\t")

WINDOW = 1000000 # 1Mb (standard window)
STATUS = c()
for(i in 1:dim(MERGED)[1])
{
  snp_gene   = MERGED$Feature_id[i]
  pheno_gene = MERGED$phenotype[i]
  if (snp_gene%in%gene_pos$gene & pheno_gene%in%gene_pos$gene)
  {
    snp_start   = gene_pos[gene_pos$gene==snp_gene,"start"]
    snp_end     = gene_pos[gene_pos$gene==snp_gene,"end"]
    pheno_start = gene_pos[gene_pos$gene==pheno_gene,"start"]
    pheno_end   = gene_pos[gene_pos$gene==pheno_gene,"end"]
    #print(paste(pheno_gene, pheno_start, pheno_end))
    status      = "trans"
    # |-1Mb/---/+1Mb|
    #       |-------–pheno-------–|
    if (((snp_start-WINDOW) < pheno_start) & ((snp_end+WINDOW) > pheno_start))
    {
      status = "cis"
    }
    #          |-1Mb/---/+1Mb|
    #       |-------–pheno-------–|
    if (((snp_start-WINDOW) > pheno_start) & ((snp_end+WINDOW) < pheno_end))
    {
      status = "cis"
    }
    #                       |-1Mb/---/+1Mb|
    #       |-------–pheno-------–|
    if (((snp_start-WINDOW) > pheno_end) & ((snp_end+WINDOW) < pheno_end))
    {
      status = "cis"
    }
  } else {
    status = "unknown"
  }
  STATUS = c(STATUS, status)
}
MERGED$snp_status = STATUS

rm(gene_pos)

#--------------------------------------------#
# 6) Collect SNP association to fitness      #
#--------------------------------------------#
SNP_FITNESS            = readRDS(paste0("./data/tribolium_eqtl/rds/",POPULATION,"_",VERSION,IMPUTATION,"_FITNESS_fitness.rds"))
names(SNP_FITNESS)[2]  = "ID"
names(SNP_FITNESS)[3]  = "pos"
names(SNP_FITNESS)[12] = "p"
SNP_FITNESS            = SNP_FITNESS[,c("ID", "beta", "p")]
names(SNP_FITNESS)     = c("ID", "snp_fitness_beta", "snp_fitness_p")
MERGED                 = merge(MERGED, SNP_FITNESS, by="ID", all.x=T)

#--------------------------------------------#
# 7) Save the resulting dataset              #
#--------------------------------------------#
saveRDS(MERGED, file=paste0("./data/tribolium_eqtl/merged/",POPULATION,"_",VERSION,IMPUTATION,"_",PHENOTYPE,"_merged.rds"))
write.table(MERGED, file=paste0("./data/tribolium_eqtl/merged/",POPULATION,"_",VERSION,IMPUTATION,"_",PHENOTYPE,"_merged.txt"), row.names=F, col.names=T, quote=F, sep="\t")

