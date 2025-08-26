#!/usr/bin/env python3
# coding: utf-8

#*******************************************************************************
# Copyright © 2021-2025 Charles Rocabert, Frédéric Guillaume
# Web: github.com/charlesrocabert/Koch-et-al-Predictability-of-Gene-Expression
#
# 4_RelocateBamFiles.py
# ---------------------
# Relocate the BAM files in a specific folder for each version.
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
import time
import argparse
import subprocess

### Parse command line arguments ###
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-path", "-repository-path", help="Repository path")
    parser.add_argument("--source-path", "-source-path", help="BAM files original source path")
    args = parser.parse_args()
    return(vars(args))

### Load the BAM map ###
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

### Write the BAM map in a specified file ###
def write_bam_map( samples, filename ):
    f      = open(filename, "w")
    header = samples[0]
    f.write(";".join(header)+"\n")
    for sample in samples[1:]:
        line = []
        for elmt in header:
            line.append(sample[elmt])
        f.write(";".join(line)+"\n")
    f.close()

### Initialize the relocation folder ###
def init_folder( path, version ):
    # Confirm with the user
    print("> You are about to erase the folder "+path+"/"+version)
    user_input = input('Confirm? [Y/N] ')
    # Input validation
    if user_input.lower() not in ('y', 'yes'):
        print("> Exit.")
        sys.exit()
    # Delete data
    if os.path.exists(path+"/"+version):
        os.system("rm -rf "+path+"/"+version)
    os.system("mkdir "+path+"/"+version)

### Relocate the files ###
def relocate_bam( samples, path, version ):
    relocated_samples = []
    relocated_samples.append(samples[0])
    for sample in samples[1:]:
        new_bam_path = path+"/"+version+"/"+sample["sample"]+".bam"
        new_bai_path = path+"/"+version+"/"+sample["sample"]+".bai"
        if not os.path.isfile(new_bam_path):
            print("> Copy BAM "+sample["sample"])
            os.system("cp "+sample["bam_path"]+" "+new_bam_path)
        sample["bam_path"] = new_bam_path
        sample["bai_path"] = new_bai_path
        relocated_samples.append(sample)
    write_bam_map(relocated_samples, "data/tribolium_bam/bam_map_ALL_"+version+"_relocated.csv")


##################
#      MAIN      #
##################

if __name__ == '__main__':

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(config["repository_path"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Relocate Tcas3.30 BAM files  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    version = "Tcas3.30"
    samples = load_bam_map("data/tribolium_bam/bam_map_ALL_"+version+".csv")
    init_folder(config["source_path"], "V1")
    relocate_bam(samples, config["source_path"], "V1")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Relocate Tcas5.2 BAM files   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    version = "Tcas5.2"
    samples = load_bam_map("data/tribolium_bam/bam_map_ALL_"+version+".csv")
    init_folder(config["source_path"], "V2")
    relocate_bam(samples, config["source_path"], "V2")

