#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 3_Gemma.py
# ----------
# Run GEMMA for a given set of phenotypes.
# (HPC SCRIPT --> array wrapper)
#
# 1) GEMMA association:
#   1.1) Import GEMMA and datasets from Allas,
#   1.2) Calculate kinship matrix,
#   1.3) Run GEMMA on the right set of phenotypes,
#   1.4) Convert the output to RDS format,
#   1.5) Export resulting files to Allas.
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
    parser.add_argument("--nb-phenotypes", "-nb-phenotypes", type=int, default=0, help="Number of phenotypes per index")
    parser.add_argument("--missing", "-missing", type=float, default=0.5, help="Missing genotypes threshold")
    parser.add_argument("--sindex", "-sindex", type=int, default=0, help="Index")
    args = parser.parse_args()
    return(vars(args))

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Import GEMMA binary ###
def import_gemma_binary( bucket ):
    os.system("a-get "+bucket+"/eqtl/gemma-0.98.5-linux-static-AMD64.gz")
    os.system("gunzip gemma-0.98.5-linux-static-AMD64.gz")
    os.system("mv gemma-0.98.5-linux-static-AMD64 gemma")
    os.system("chmod ug+x gemma")
    assert_creation("gemma")

### Import the dataset ###
def import_data( bucket, dataset ):
    os.system("a-get "+bucket+"/eqtl/"+dataset+".bed")
    os.system("a-get "+bucket+"/eqtl/"+dataset+".bim")
    os.system("a-get "+bucket+"/eqtl/"+dataset+".fam")
    os.system("a-get "+bucket+"/eqtl/"+dataset+".pheno")
    assert_creation(dataset+".bed")
    assert_creation(dataset+".bim")
    assert_creation(dataset+".fam")
    assert_creation(dataset+".pheno")

### Load the list of phenotypes ###
def load_phenotypes( dataset, nb_phenotypes ):
    phenotypes  = {1:[]}
    sindex      = 1
    count       = 0
    pheno_index = 1
    f = open(dataset+".pheno", "r")
    l = f.readline()
    while l:
        pheno_name = l.strip("\n")
        if count < nb_phenotypes:
            phenotypes[sindex].append([pheno_name, pheno_index])
            count       += 1
            pheno_index += 1
        else:
            count   = 1
            sindex += 1
            phenotypes[sindex] = [[pheno_name, pheno_index]]
            pheno_index += 1
        l = f.readline()
    f.close()
    return phenotypes

### Compute the kinship matrix ###
def compute_kinship_matrix( dataset ):
    cmdline_elmts = []
    cmdline_elmts.append("./gemma")
    cmdline_elmts.append("-bfile "+dataset)
    cmdline_elmts.append("-gk")
    cmdline_elmts.append("-o "+dataset)
    cmdline = " ".join(cmdline_elmts)
    os.system(cmdline)
    assert_creation("output/"+dataset+".cXX.txt")

### Run GEMMA association for a given phenotype ###
def run_gemma( dataset, missing, pheno_name, pheno_index ):
    cmdline_elmts = []
    cmdline_elmts.append("./gemma")
    cmdline_elmts.append("-bfile "+dataset)
    cmdline_elmts.append("-k ./output/"+dataset+".cXX.txt")
    cmdline_elmts.append("-lmm")
    cmdline_elmts.append("-n "+str(pheno_index))
    cmdline_elmts.append("-miss "+str(missing))
    cmdline_elmts.append("-notsnp")
    cmdline_elmts.append("-o "+dataset+"_"+pheno_name)
    cmdline = " ".join(cmdline_elmts)
    print(cmdline)
    os.system(cmdline)
    assert_creation("output/"+dataset+"_"+pheno_name+".assoc.txt")
    assert_creation("output/"+dataset+"_"+pheno_name+".log.txt")
    os.system("gzip -c output/"+dataset+"_"+pheno_name+".assoc.txt > output/"+dataset+"_"+pheno_name+".assoc.txt.gz")

### Create the conversion R-script ###
def convert_to_rds( bucket, dataset, pheno_name, local_path ):
    gemma_file = "output/"+dataset+"_"+pheno_name+".assoc.txt.gz"
    rds_file   = "output/"+dataset+"_"+pheno_name+".rds"
    f          = open("ConvertToRDS.R", "w")
    f.write('print("'+local_path+'")\n')
    f.write('setwd("'+local_path+'")\n')
    f.write('f = "'+gemma_file+'"\n')
    f.write('print(f)\n')
    f.write('d = read.table(gzfile(f), h=T, sep="\t", check.names=F)\n')
    f.write('print(head(d))\n')
    f.write('saveRDS(d, file="'+rds_file+'")\n')
    f.close()
    os.system("singularity_wrapper exec Rscript ConvertToRDS.R")
    assert_creation(rds_file)
    os.system("rm "+gemma_file)

### Export resulting files ###
def export_files( bucket, dataset, pheno_name ):
    rds_file   = dataset+"_"+pheno_name+".rds"
    os.system("cp output/"+rds_file+" /scratch/project_XXXXXXX/Tribolium_castaneum_eQTL/rds/"+rds_file)


##################
#      MAIN      #
##################

if __name__ == '__main__':
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config     = parse_arguments()
    local_path = os.environ['DATADIR']
    os.chdir(os.environ['DATADIR'])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Import GEMMA binary          #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Import GEMMA binary")
    import_gemma_binary(config["bucket"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Import the dataset           #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Import the dataset")
    import_data(config["bucket"], config["dataset"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Load the list of phenotypes  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load the list of phenotypes")
    phenotypes = load_phenotypes(config["dataset"], config["nb_phenotypes"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 5) Compute the kinship matrix   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Compute the kinship matrix")
    compute_kinship_matrix(config["dataset"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 6) Run GEMMA for each phenotype #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Run GEMMA for each phenotype")
    for item in phenotypes[config["sindex"]]:
        pheno_name  = item[0]
        pheno_index = item[1]
        print("   > Running GEMMA for phenotype "+pheno_name+" ("+str(pheno_index)+")")
        run_gemma(config["dataset"], config["missing"], pheno_name, pheno_index)
        convert_to_rds(config["bucket"], config["dataset"], pheno_name, local_path)
        export_files(config["bucket"], config["dataset"], pheno_name)

    print(">> Done.")

