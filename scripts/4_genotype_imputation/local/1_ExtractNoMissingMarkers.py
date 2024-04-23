#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_ExtractNoMissingMarkers.py
# ----------------------------
# Extract all markers with a call rate of 100%, and save them in an ad hoc
# binary data structure.
# (LOCAL SCRIPT)
#***************************************************************************

import os
import sys
import dill
import copy
import argparse

### Parse command line arguments ###
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-path", "-repository-path", help="Repository path")
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--suffix", "-suffix", help="Suffix")
    args = parser.parse_args()
    return(vars(args))

### Load the VCF file as a matrix ###
def load_VCF( filename ):
    genotypes = {"header":"", "missing":[]}
    f         = open(filename, "r", encoding="ISO-8859-1")
    ### Skip header lines ###
    l = f.readline()
    while "#CHROM" not in l:
        genotypes["header"] += l
        l = f.readline()
    ### Parse the header ###
    global_params        = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    header               = l.strip("\n").split("\t")
    genotypes["header"] += l
    l                    = f.readline()
    counter              = 1
    while l:
        if counter%10000 == 0:
            print("   > "+str(counter)+" variants parsed")
        l  = l.strip("\n").split("\t")
        ID = l[2]
        if ID == ".":
            ID = l[0]+"-"+l[1]
        assert ID not in genotypes.keys(), ID
        genotypes[ID] = {"global_params":{}}
        missing       = 0.0
        N             = 0.0
        for i in range(len(header)):
            if header[i] not in global_params:
                geno = l[i].split(":")[0]
                assert geno in ["./.", "0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1"], geno
                genotypes[ID][header[i]] = geno
                N += 1.0
                if geno == "./.":
                    missing += 1.0
            else:
                genotypes[ID]["global_params"][header[i]] = l[i]
        genotypes[ID]["missing"] = missing/N
        genotypes["missing"].append(missing)
        l = f.readline()
        counter += 1
    f.close()
    return genotypes

### Write in a file the probability distribution of missing genotypes ###
def write_probability_distribution( filename, genotypes ):
    f = open(filename, "w")
    for missing in genotypes["missing"]:
        f.write(str(missing)+"\n")
    f.close()

### Select all variants with no missing genotypes ###
def keep_no_missing_variants( genotypes ):
    no_missing            = {}
    no_missing["header"]  = genotypes["header"]
    no_missing["missing"] = genotypes["missing"][:]
    counter               = 0
    for ID in genotypes.keys():
        if ID not in ["header", "missing"]:
            if genotypes[ID]["missing"] == 0.0:
                no_missing[ID] = copy.deepcopy(genotypes[ID])
                counter += 1
    print("   > "+str(counter)+" SNPs selected.")
    return no_missing


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
    # 2) Load the dataset             #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load the VCF file")
    filename  = "./data/tribolium_snp/Tribolium_castaneum_"+config["population"]+"_"+config["version"]+"_"+config["suffix"]+".vcf"
    genotypes = load_VCF(filename)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Save the distribution of     #
    #    missing genotypes            #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Save the probability distribution of missing genotypes")
    write_probability_distribution("./data/tribolium_snp/imputation_tests/missing_genotypes_probability_distribution.csv", genotypes)
    print(">> Select SNPs with 100% call rate")
    no_missing = keep_no_missing_variants(genotypes)
    del genotypes

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Dump the result in a binary  #
    #    file                         #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ofile = open("./data/tribolium_snp/imputation_tests/no_missing_genotypes", "wb")
    dill.dump(no_missing, ofile)
    ofile.close()
    del no_missing

