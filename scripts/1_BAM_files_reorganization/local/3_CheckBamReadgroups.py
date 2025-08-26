#!/usr/bin/env python3
# coding: utf-8

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# 3_CheckBamReadgroups.py
# -----------------------
# Check the absence of the readgroup in every Tcas3.30 (2016) and Tcas5.2
# (2017) BAM files.
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

import os
import sys
import csv
import argparse
import subprocess

### Parse command line arguments ###
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-path", "-repository-path", help="Repository path")
    args = parser.parse_args()
    return(vars(args))

### Load the bam map ###
def load_bam_map( filename ):
    samples   = []
    file      = open(filename, "r")
    csvreader = csv.reader(file, delimiter=";")
    header    = next(csvreader)
    samples.append(header)
    for row in csvreader:
        sample = {}
        for i in range(len(header)):
            sample[header[i]] = row[i]
        samples.append(sample)
    file.close()
    return samples

### Check the existence of a read group ###
def check_readgroup( samples ):
    for sample in samples[1:]:
        bam_path = sample["bam_path"]
        proc     = subprocess.Popen(["samtools view -H "+bam_path+" | grep ^\@RG"], stdout=subprocess.PIPE, shell=True)
        output   = proc.stdout.read().decode('utf8')
        if len(output.strip("\n")) > 0:
            print("> sample "+sample["sample"]+" "+sample["genome_ref"]+" "+sample["annotation"])
            print("  >>> Read group is present.")


##################
#      MAIN      #
##################

if __name__ == '__main__':

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(config["repository_path"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Check Tcas3.30 BAM read groups #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print("> Evaluate Tcas3.30 (2016) samples")
    tcas3_30_samples = load_bam_map("data/tribolium_bam/bam_map_ALL_Tcas3.30.csv")
    check_readgroup(tcas3_30_samples)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Check Tcas5.2 BAM read groups  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print("> Evaluate Tcas5.2 (2017) samples")
    tcas5_2_samples = load_bam_map("data/tribolium_bam/bam_map_ALL_Tcas5.2.csv")
    check_readgroup(tcas5_2_samples)

