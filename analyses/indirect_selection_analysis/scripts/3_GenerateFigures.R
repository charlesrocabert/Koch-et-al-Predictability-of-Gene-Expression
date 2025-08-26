#!/usr/bin/env Rscript

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# 3_GenerateFigures.R
# -------------------
# Generate all figures associated to the manuscript.
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
library("cowplot")
library("ggpubr")
library("ggbreak")
library("patchwork")
library("RColorBrewer")
library("rstudioapi")
library("latex2exp")

### Plot gene parallelism against total selection ###
gene_parallelism_vs_selection <- function( gene_dataset )
{
  ############
  comparisons = list(c("No shift", "Single shift"),
                     c("No shift", "Partial/full parallelism"),
                     c("Single shift", "Partial/full parallelism"))
  ############
  count_summary = gene_dataset %>% group_by(Parallel_category) %>% tally()
  ############
  p = ggplot(gene_dataset, aes(Parallel_category, log_total_selection, fill=Parallel_category)) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary$n), x=c(1,2,3), y=rep(3.5,3), size=3, vjust=2) +
    ylim(-6, 6.2) +
    xlab("") +
    ylab(TeX("Genetic selection - $|\\log_{10}|$")) +
    ggtitle("Gene AFC parallelism and selection") +
    labs(fill = "Parallelism:") +
    theme_classic()
  ############
  return(p)
}

### Plot eQTL-SNPs parallelism against total selection ###
eQTL_parallelism_vs_selection <- function( eQTL_dataset )
{
  ############
  comparisons = list(c("No shift", "Single shift"),
                     c("No shift", "Partial/full parallelism"),
                     c("Single shift", "Partial/full parallelism"))
  ############
  count_summary = eQTL_dataset[!duplicated(eQTL_dataset$ID),] %>% group_by(Parallel_category) %>% tally()
  ############
  p = ggplot(eQTL_dataset[!duplicated(eQTL_dataset$ID),], aes(Parallel_category, log_total_selection, fill=Parallel_category)) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary$n), x=c(1,2,3), y=rep(1.6,3), size=3, vjust=2) +
    ylim(-6, 4) +
    xlab("") +
    ylab("Weighted total selection\n(absolute log-scale)") +
    ggtitle("eQTLs") +
    labs(fill = "Parallelism:") +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  return(p)
}

### Plot eQTL carrier categories against total selection ###
eQTL_carrier_vs_selection <- function( gene_dataset )
{
  ############
  comparisons = list(c("Non-eQTL carrier", "Single-phenotype eQTL carrier"),
                     c("Non-eQTL carrier", "Pleiotropic eQTL carrier"),
                     c("Single-phenotype eQTL carrier", "Pleiotropic eQTL carrier"))
  ############
  count_summary = gene_dataset %>% group_by(Pleiotropy_category) %>% tally()
  ############
  p = ggplot(gene_dataset, aes(Pleiotropy_category, log_total_selection, fill=Pleiotropy_category)) +
    geom_boxplot() +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary$n), x=c(1,2,3), y=rep(3.5,3), size=3, vjust=2) +
    scale_fill_brewer(palette="PiYG", labels=c("Non-eQTL\ncarrier", "Single-phenotype\neQTL carrier", "Pleiotropic\neQTL carrier")) +
    scale_x_discrete(labels=c("Non-eQTL\ncarrier", "Single-phenotype\neQTL carrier", "Pleiotropic\neQTL carrier")) +
    xlab("") +
    ylab(TeX("Genetic selection - $|\\log_{10}|$")) +
    ggtitle("eQTL carrier categories and selection") +
    labs(fill = "") +
    theme_classic()
  ############
  return(p)
}

### Plot eQTL phenotype categories against total selection ###
eQTL_phenotype_vs_selection <- function( gene_dataset )
{
  ############
  comparisons = list(c("Non-eQTL phenotype", "Low connectivity"),
                     c("Non-eQTL phenotype", "High connectivity"),
                     c("Low connectivity", "High connectivity"))
  ############
  count_summary = gene_dataset %>% group_by(Connectivity_category) %>% tally()
  ############
  p = ggplot(gene_dataset, aes(Connectivity_category, log_total_selection, fill=Connectivity_category)) +
    geom_boxplot() +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary$n), x=c(1,2,3), y=rep(3,3), size=3, vjust=2) +
    scale_fill_brewer(palette="PuOr", labels=c("Non-eQTL\nphenotype", "Low\nconnectivity (<5)", "High\nconnectivity (>4)")) +
    xlab("") +
    ylab("Total selection on gene expression\n(absolute log-scale)") +
    ggtitle("eQTL phenotype categories and selection") +
    labs(fill = "") +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  return(p)
}

### Plot hub genes against total selection ###
hub_gene_figures <- function( gene_dataset )
{
  ############
  comparisons = list(c("Other","Hub gene"))
  ############
  count_summary = gene_dataset %>% group_by(hub_gene_category) %>% tally()
  ############
  p = ggplot(gene_dataset, aes(hub_gene_category, log_total_selection, fill=hub_gene_category)) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary$n), x=c(1,2), y=rep(3,2), size=3, vjust=2) +
    xlab("") +
    ylab("Total selection on gene expression\n(absolute log-scale)") +
    ggtitle("Hub genes selection compared to other genes") +
    labs(fill = "") +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  return(p)
}

### Distribution of eQTL pleiotropy and eQTL phenotype connectivity ###
eQTL_histograms <- function( eQTL_carrier_dataset, eQTL_phenotype_dataset )
{
  ############
  p1 = ggplot(eQTL_carrier_dataset, aes(x=Pleiotropy)) +
    geom_histogram(fill="cornflowerblue", alpha=0.5, color="black") +
    scale_y_cut(breaks=c(25), which=c(1), scales=c(0.5)) +
    xlab("Pleiotropy of eQTL carrier genes") +
    ylab("Count") +
    ggtitle("Distribution of eQTL carrier genes' pleiotropy") +
    theme_classic()
  ############
  p2 = ggplot(eQTL_phenotype_dataset, aes(x=Connectivity)) +
    geom_histogram(fill="cornflowerblue", alpha=0.5, color="black") +
    scale_y_cut(breaks=c(100), which=c(1), scales=c(0.5)) +
    xlab("Phenotypes connectivity") +
    ylab("Count") +
    ggtitle("Distribution of phenotypes' connectivity") +
    theme_classic()
  ############
  #p1 = plot_grid(p1, labels=c("A"))
  #p2 = plot_grid(p2, labels=c("B"))
  #p  = plot_grid(p1, p2)
  return(list("p1"=p1, "p2"=p2))
  #return(p)
}

### Distribution of Jaccard similarity index between lines AFCs ###
line_shift_similarity_figure <- function( SNP_dataset )
{
  LINES = c("L1_shift", "L2_shift", "L3_shift", "L5_shift", "L6_shift", "Mx1_shift", "Mx2_shift")
  D     = SNP_dataset[,LINES]
  S     = matrix(rep(0,7*7), ncol=7)
  N     = nrow(D)
  Z     = c()
  for(i in seq(1,6))
  {
    for(j in seq(i+1,7))
    {
      d = D[,c(i,j)]
      x = rowSums(d)
      d = d[which(x!=0),]
      n = nrow(d)
      x = rowSums(d)
      z = sum(x==2)
      Z = c(Z, z/n*100)
    }
  }
  Z = data.frame(Z)
  p = ggplot(Z, aes(x=0, y=Z)) +
    geom_violin(fill=rgb(210/255, 180/255, 112/255)) +
    geom_boxplot(width=0.1, color="black", fill="white", alpha=1) +
    theme_classic() +
    ggtitle("Distribution of frequency shift similarities between lines") +
    xlab("") +
    ylab("Similarity\n(% of SNPs shifting in common between two lines)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  print(paste0("> Mean similarity = ",mean(Z$Z)))
  print(paste0("> Similarity s.e. = ",sd(Z$Z)))
  return(p)
}

### Gene expression VS. number of SNPs per gene ###
expression_vs_nbSNPs_figure <- function( gene_dataset )
{
  dl  = filter(gene_dataset, avg_read_count > 0)
  reg = lm(log10(nb_SNPs/gene_length) ~ log10(avg_read_count), data=dl)
  res = summary(reg)
  p   = ggplot(dl, aes(x=avg_read_count, y=nb_SNPs/gene_length)) +
    geom_point() +
    geom_smooth() +
    xlab("Average read count") +
    ylab("Number of SNPs\n(normalized by gene length)") +
    ggtitle("Number of SNPs per gene VS. average read count") +
    scale_x_log10() + scale_y_log10() +
    annotate("text", x=1e+5, y=1e-5, label=paste0("R² = ", round(res$r.squared, 3)), size=5) +
    theme_classic()
  return(p)
}

### LD decay with distance ###
LD_decay_figure <- function( gene_dataset, LD_decay_dataset, window )
{
  LD_decay_dataset$distc = cut(LD_decay_dataset$dist,breaks=seq(from=min(LD_decay_dataset$dist)-1,to=max(LD_decay_dataset$dist)+1,by=window))
  Davg  = LD_decay_dataset %>% group_by(distc) %>% summarise(mean=mean(R2),median=median(R2))
  Davg  = Davg %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                            end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                            mid=start+((end-start)/2))
  Davg  = filter(Davg, !is.na(mid))
  Gmean = mean(log10(gene_dataset$gene_length))
  Gsdm  = mean(log10(gene_dataset$gene_length))-sd(log10(gene_dataset$gene_length))
  Gsdp  = mean(log10(gene_dataset$gene_length))+sd(log10(gene_dataset$gene_length))
  p1 = ggplot(Davg, aes(x=mid/1000000, y=mean)) +
    geom_point() +
    geom_line() +
    geom_vline(aes(xintercept=10^Gmean/1000000, color="Average gene length"), linetype="dashed") +
    geom_vline(aes(xintercept=10^Gsdm/1000000, color="Average +/- SD"), linetype="twodash") +
    geom_vline(aes(xintercept=10^Gsdp/1000000, color="Average +/- SD"), linetype="twodash") +
    xlab("Distance (Mb)") +
    ylab("LD (r²)") +
    ggtitle("LD decay with distance (in Mb)") +
    scale_x_log10() + #scale_y_log10() +
    scale_color_discrete(name = "Legend") +
    theme_classic() +
    theme(legend.position="bottom")
  p2 = ggplot(gene_dataset, aes(x=(gene_length/1000000))) +
    geom_histogram(binwidth=0.1) +
    xlab("Gene length (Mb)") +
    ylab("Number of genes") +
    ggtitle("Gene length distribution") +
    scale_x_log10() +
    theme_classic()
  p = plot_grid(p1, p2, ncol=2, labels="AUTO")
  return(p)
}

### Plot the number of SNPs in hub genes vs. other genes ###
hub_nb_SNPs_figure <- function( gene_dataset )
{
  gene_dataset$hub_category = as.character(gene_dataset$hub_gene)
  gene_dataset$hub_category[gene_dataset$hub_category=="0"] = "Other"
  gene_dataset$hub_category[gene_dataset$hub_category=="1"] = "Hub gene"
  gene_dataset$hub_category = factor(gene_dataset$hub_category, levels=c("Other","Hub gene"))
  comparisons = list(c("Other","Hub gene"))
  count_summary = gene_dataset %>% group_by(hub_category) %>% tally()
  p1 = ggplot(gene_dataset, aes(as.factor(hub_category), log10(nb_SNPs), fill=as.factor(hub_category))) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary$n), x=c(1,2), y=rep(0.2,2), size=4, vjust=2) +
    #ylim(-6, 6.2) +
    xlab("") +
    ylab(TeX("Number of SNPs - $|\\log_{10}|$")) +
    ggtitle("Number of SNPs in hub or other genes") +
    labs(fill = "Gene category:") +
    theme_classic()
  p2 = ggplot(gene_dataset, aes(as.factor(hub_category), log10(gene_length), fill=as.factor(hub_category))) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary$n), x=c(1,2), y=rep(0.2,2), size=4, vjust=2) +
    #ylim(-6, 6.2) +
    xlab("") +
    ylab(TeX("Gene length - $|\\log_{10}|$")) +
    ggtitle("Number of SNPs in hub or other genes") +
    labs(fill = "Gene category:") +
    theme_classic()
  p3 = ggplot(gene_dataset, aes(as.factor(hub_category), log10(nb_SNPs/gene_length), fill=as.factor(hub_category))) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary$n), x=c(1,2), y=rep(0.5,2), size=4, vjust=2) +
    #ylim(-6, 6.2) +
    xlab("") +
    ylab(TeX("Number of SNPs (normalized by gene length) - $|\\log_{10}|$")) +
    ggtitle("Number of SNPs in hub or other genes") +
    labs(fill = "Gene category:") +
    theme_classic()
  return(list("p1"=p1, "p2"=p2, "p3"=p3))
}


##################
#      MAIN      #
##################

directory = dirname(getActiveDocumentContext()$path)
setwd(directory)
setwd("../../..")

#------------------------------------------------#
# 1) Load the data                               #
#------------------------------------------------#
SNP_dataset            = readRDS("./analyses/indirect_selection_analysis/data/SNP_dataset.rds")
gene_dataset           = readRDS("./analyses/indirect_selection_analysis/data/gene_dataset.rds")
eQTL_dataset           = readRDS("./analyses/indirect_selection_analysis/data/eQTL_dataset.rds")
eQTL_carrier_dataset   = readRDS("./analyses/indirect_selection_analysis/data/eQTL_carrier_dataset.rds")
eQTL_phenotype_dataset = readRDS("./analyses/indirect_selection_analysis/data/eQTL_phenotype_dataset.rds")
LD_decay_dataset       = readRDS("./analyses/indirect_selection_analysis/data/LD_decay_dataset.rds")

#------------------------------------------------#
# 2) Build the main figure                       #
#------------------------------------------------#
p1 = gene_parallelism_vs_selection(gene_dataset) + theme(legend.position="none")
p2 = eQTL_carrier_vs_selection(gene_dataset) + theme(legend.position="none")
p  = plot_grid(p1, p2, labels="AUTO")
ggsave("./analyses/indirect_selection_analysis/plots/genetic_markers_total_selection.pdf", p, width=9, height=5, units="in")

#------------------------------------------------#
# 3) Make pleiotropy and connectivity histograms #
#------------------------------------------------#
figs = eQTL_histograms(eQTL_carrier_dataset, eQTL_phenotype_dataset)
ggsave("./analyses/indirect_selection_analysis/plots/pleiotropy_distribution.pdf", figs[["p1"]], width=5, height=4, units="in")
ggsave("./analyses/indirect_selection_analysis/plots/connectivity_distribution.pdf", figs[["p2"]], width=5, height=4, units="in")

#------------------------------------------------#
# 4) Make hub gene figures                       #
#------------------------------------------------#
p = hub_gene_figures(gene_dataset)
ggsave("./analyses/indirect_selection_analysis/plots/hub_genes_total_selection.pdf", p, width=5, height=4, units="in")

#------------------------------------------------#
# 5) Make the expression VS. SNP number figure   #
#------------------------------------------------#
p = expression_vs_nbSNPs_figure(gene_dataset)
ggsave("./analyses/indirect_selection_analysis/plots/expression_vs_nb_SNPs.pdf", p, width=8, height=4, units="in")

#------------------------------------------------#
# 6) Make the LD decay figure                    #
#------------------------------------------------#
p = LD_decay_figure(gene_dataset, LD_decay_dataset, 1500)
ggsave("./analyses/indirect_selection_analysis/plots/LD_decay.pdf", p, width=9, height=4, units="in")

#------------------------------------------------#
# 7) Make nb SNPs in hub genes figures           #
#------------------------------------------------#
figs = hub_nb_SNPs_figure(gene_dataset)
ggsave("./analyses/indirect_selection_analysis/plots/hub_genes_nb_SNPs.pdf", figs[["p3"]], width=5, height=4, units="in")

# x = gene_dataset$nb_SNPs
# y = c()
# for(elmt in x)
# {
#   parallel = c()
#   for(i in seq(1, elmt))
#   {
#     res      = rbinom(7, 1, 0.025)
#     parallel = c(parallel, sum(res))
#   }
#   y   = c(y, max(parallel))
# }
# sim_data        = data.frame(x, y)
# names(sim_data) = c("nb_SNPs", "Parallelism")
# sim_data$source = rep("Simulation", nrow(sim_data))
# D = data.frame(gene_dataset$nb_SNPs, gene_dataset$Parallelism)
# names(D) = c("nb_SNPs", "Parallelism")
# D$source = rep("Experimental", nrow(D))
# D = rbind(D, sim_data)
# ggplot(D, aes(as.factor(Parallelism), nb_SNPs, fill=source)) +
#   geom_boxplot() +
#   scale_y_log10()
# 
# X = as.vector(table(gene_dataset$Parallelism))
# Y = as.vector(table(sim_data$Parallelism))
# 
# plot(X, type="l", ylim=range(c(X,Y)))
# lines(Y)

