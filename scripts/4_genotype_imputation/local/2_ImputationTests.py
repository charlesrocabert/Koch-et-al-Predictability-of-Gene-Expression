#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 2_ImputationTests.py
# --------------------
# Test imputation capabilities of Beagle with benchmark datasets.
# (LOCAL SCRIPT)
#***************************************************************************

import os
import sys
import csv
import copy
import time
import dill
import random
import argparse

### Parse command line arguments ###
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-path", "-repository-path", help="Repository path")
    parser.add_argument("--beagle", "-beagle", help="Beagle-5.4 path")
    parser.add_argument("--rep", "-rep", type=int, help="Number of repetitions")
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

### Generate a toy dataset with random missing genotypes ###
def generate_toy_data( genotypes ):
    toy_data            = {}
    toy_data["header"]  = genotypes["header"]
    toy_data["missing"] = []
    for ID in genotypes.keys():
        if ID not in ["header", "missing"]:
            toy_data[ID]            = copy.deepcopy(genotypes[ID])
            N_missing               = random.choice(genotypes["missing"])
            N                       = float(len(toy_data[ID])-2)
            toy_data[ID]["missing"] = N_missing/N
            toy_data["missing"].append(N_missing)
            pos_vec  = range(int(N))
            pos_miss = random.choices(pos_vec, k=int(N_missing))
            pos      = 0
            for key in toy_data[ID].keys():
                if key not in ["missing", "global_params"]:
                    if pos in pos_miss:
                        toy_data[ID][key] = "./."
                    pos += 1
    return toy_data

### Write the toy data into a VCF ###
def write_toy_VCF( toy_data ):
    global_params = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    f             = open("./data/tribolium_snp/imputation_tests/toy_data.vcf", "w", encoding="ISO-8859-1")
    f.write(toy_data["header"])
    for ID in toy_data.keys():
        if ID not in ["header", "missing"]:
            line = ""
            for param in global_params:
                line += toy_data[ID]["global_params"][param]+"\t"
            for key in toy_data[ID].keys():
                if key not in ["missing", "global_params"]:
                    line += toy_data[ID][key]+"\t"
            f.write(line.strip("\t")+"\n")
    f.close()

### Run Beagle for imputation ###
def run_beagle_imputation( beagle_path ):
    cmdline = "java -jar "+beagle_path+" gt=./data/tribolium_snp/imputation_tests/toy_data.vcf out=./data/tribolium_snp/imputation_tests/imputed_data"
    os.system(cmdline+" > /dev/null")
    os.system("gunzip ./data/tribolium_snp/imputation_tests/imputed_data.vcf.gz")

### Load the imputed VCF file ###
def compute_distance( complete_data, missing_data, imputed ):
    SUCCESS = 0.0
    COUNT   = 0.0
    for ID in complete_data.keys():
        if ID not in ["header", "missing"]:
            assert ID in missing_data.keys(), ID
            assert ID in imputed.keys(), ID
            for key in complete_data[ID].keys():
                if key not in ["missing", "global_params"]:
                    assert key in missing_data[ID].keys(), key
                    assert key in imputed[ID].keys(), key
                    geno1 = complete_data[ID][key]
                    geno2 = missing_data[ID][key]
                    geno3 = imputed[ID][key]
                    if geno2 == "./.":
                        COUNT += 1.0
                        if "/" in geno1:
                            if geno1 == "0/0" and geno3 == "0|0":
                                SUCCESS += 1.0
                            elif geno1 in ["0/1", "1/0"] and geno3 in ["0|1", "1|0"]:
                                SUCCESS += 1.0
                            elif geno1 == "1/1" and geno3 == "1|1":
                                SUCCESS += 1.0
                            #else:
                            #    print(">>> Well well ... "+geno1+" "+geno2+" "+geno3)
                        elif "|" in geno1:
                            if geno1 == "0|0" and geno3 == "0|0":
                                SUCCESS += 1.0
                            elif geno1 in ["0|1", "1|0"] and geno3 in ["0|1", "1|0"]:
                                SUCCESS += 1.0
                            elif geno1 == "1|1" and geno3 == "1|1":
                                SUCCESS += 1.0
                            #else:
                            #    print(">>> Well well ... "+geno1+" "+geno2+" "+geno3)
                        #else:
                        #    print(">>> Well well ... "+geno1+" "+geno2+" "+geno3)
    return SUCCESS, COUNT


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
    print(">> Load the 100% call rate dataset")
    ifile      = open("./data/tribolium_snp/imputation_tests/no_missing_genotypes", "rb")
    no_missing = dill.load(ifile)
    ifile.close()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Run imputation on toy data   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Run the tests")
    random.seed(time.time())
    f = open("./data/tribolium_snp/imputation_tests/imputation_success_rate.csv", "w")
    f.write("success_rate;nb_success;total\n")
    f.flush()
    for rep in range(1,config["rep"]+1):
        print("   > Repetition "+str(rep)+"...")
        print("     • Generate toy data")
        toy_data = generate_toy_data(no_missing)
        print("     • Write VCF file")
        write_toy_VCF(toy_data)
        print("     • Run genotype imputation")
        run_beagle_imputation(config["beagle"])
        imputed        = load_VCF("./data/tribolium_snp/imputation_tests/imputed_data.vcf")
        print("     • Compute success rate")
        success, count = compute_distance(no_missing, toy_data, imputed)
        print("     • Score = "+str(success/count))
        f.write(str(success/count)+";"+str(success)+";"+str(count)+"\n")
        f.flush()
        os.system("rm ./data/tribolium_snp/imputation_tests/imputed_data.vcf")
        os.system("rm ./data/tribolium_snp/imputation_tests/toy_data.vcf")
    f.close()

    print(">> Done.")

