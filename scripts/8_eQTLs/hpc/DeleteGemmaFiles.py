#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# DeleteGemmaFiles.py
# -------------------
# Delete GEMMA output files.
# (HPC SCRIPT -> run wrapper)
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
    parser.add_argument("--bucket", "-bucket", help="Allas bucket name")
    parser.add_argument("--dataset", "-dataset", help="Dataset name")
    args = parser.parse_args()
    return(vars(args))

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Import the phenotypes ###
def import_phenotypes( bucket, dataset ):
    os.system("a-get "+bucket+"/eqtl/"+dataset+".pheno")
    assert_creation(dataset+".pheno")

### Load the list of phenotypes ###
def load_phenotypes( dataset ):
    phenotypes = []
    f = open(dataset+".pheno", "r")
    l = f.readline()
    while l:
        pheno_name = l.strip("\n")
        phenotypes.append(pheno_name)
        l = f.readline()
    f.close()
    return phenotypes

### Delete GEMMA files for a given phenotype ###
def delete_GEMMA_files( bucket, dataset, pheno_name ):
    os.system("a-delete -f "+bucket+"/eqtl/rds/"+dataset+"_"+pheno_name+".rds")


##################
#      MAIN      #
##################

if __name__ == '__main__':
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments           #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Import and load the list of phenotypes #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Import and load the list of phenotypes")
    import_phenotypes(config["bucket"], config["dataset"])
    phenotypes = load_phenotypes(config["dataset"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Delete output files for each phenotype #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Delete output files for each phenotype")
    for i in range(len(phenotypes)):
        if i%100 == 0:
            print("   > "+str(i+1)+" phenotypes scanned")
        pheno_name = phenotypes[i]
        delete_GEMMA_files(config["bucket"], config["dataset"], pheno_name)
        time.sleep(0.1)
    os.system("rm "+config["dataset"]+".pheno")

    print(">> Done.")

