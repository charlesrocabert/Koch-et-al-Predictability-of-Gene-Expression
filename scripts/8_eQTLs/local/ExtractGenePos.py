#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# ExtractGenePos.py
# -----------------
# Extract gene positions from the GFF file.
# (LOCAL SCRIPT)
#***************************************************************************

import os
import sys
import csv
import time
import argparse
import subprocess

### Parse command line arguments ###
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-path", "-repository-path", help="Path to the repository")
    parser.add_argument("--version", "-version", help="Reference genome version")
    args = parser.parse_args()
    return(vars(args))

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Parse the GFF file ###
def parse_gff_file( version ):
    gene_pos = {}
    f        = open("./data/tribolium_genome/Tribolium_castaneum_"+version+"/Tribolium_castaneum_"+version+".gff3", "r")
    l        = f.readline()
    while l:
        if not l.startswith("#"):
            l       = l.strip("\n").split("\t")
            chr     = l[0]
            feature = l[2]
            if chr != "Unknown" and feature in ["gene", "pseudogene"]:
                gene_name = l[8].strip("ID=gene:").split(";")[0]
                assert gene_name not in gene_pos.keys()
                gene_pos[gene_name] = [chr, int(l[3]), int(l[4])]
        l = f.readline()
    f.close()
    return gene_pos


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
    # 2) Load the list of phenotypes  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load the list of gene positions")
    gene_pos = parse_gff_file(config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Save the resulting dataset   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Save the resulting dataset")
    f = open("./data/tribolium_eqtl/gene_pos_"+config["version"]+".csv", "w")
    f.write("gene\tchr\tstart\tend\n")
    for item in gene_pos.items():
        f.write(item[0]+"\t"+item[1][0]+"\t"+str(item[1][1])+"\t"+str(item[1][2])+"\n")
    f.close()

    print(">> Done.")

