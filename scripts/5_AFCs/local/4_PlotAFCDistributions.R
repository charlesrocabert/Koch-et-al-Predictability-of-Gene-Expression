#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 4_PlotAFCDistributions.R
# ------------------------
# Plot AFC distributions.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")

##################
#      MAIN      #
##################

#--------------------------------#
# 1) Read command line arguments #
#--------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
setwd(REPOSITORY_PATH)

#--------------------------------#
# 2) Load the AFC datasets       #
#--------------------------------#
AFC_CT = readRDS("./data/tribolium_afc/AFC_CT.rds")
AFC_HD = readRDS("./data/tribolium_afc/AFC_HD.rds")

#--------------------------------#
# 3) Build some plots            #
#--------------------------------#
ggplot() +
  geom_density(data=filter(AFC_CT, AFC!=0.0), aes(abs(AFC), color="CT")) +
  geom_density(data=filter(AFC_HD, AFC!=0.0), aes(abs(AFC), color="HD")) +
  scale_y_log10()

ggplot() +
  geom_density(data=AFC_CT, aes(abs(AFC), color="CT")) +
  geom_density(data=AFC_HD, aes(abs(AFC), color="HD")) +
  scale_y_log10()
