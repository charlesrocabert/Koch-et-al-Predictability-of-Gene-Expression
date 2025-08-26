#!/usr/bin/env Rscript

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# 4_PlotAFCDistributions.R
# ------------------------
# Plot AFC distributions.
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
