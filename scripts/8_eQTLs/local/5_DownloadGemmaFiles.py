#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 5_DownloadGemmaFiles.py
# -----------------------
# Download Gemma association files from a distant server.
# (LOCAL SCRIPT)
#***************************************************************************

import os
import sys
import csv
import time
import argparse
import subprocess
import swiftclient
import paramiko

### Parse command line arguments ###
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-path", "-repository-path", help="Path to the repository")
    parser.add_argument("--user", "-user", help="Puhti user name")
    parser.add_argument("--password", "-password", help="Puhti password")
    parser.add_argument("--dataset", "-dataset", help="Dataset name")
    args = parser.parse_args()
    return(vars(args))

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Load the list of phenotypes ###
def load_phenotypes( dataset ):
    phenotypes = []
    f = open("./data/tribolium_eqtl/gemma/"+dataset+".pheno", "r")
    l = f.readline()
    while l:
        pheno_name = l.strip("\n")
        phenotypes.append(pheno_name)
        l = f.readline()
    f.close()
    return phenotypes

### Connect to Puhti via SSH ###
def connect_to_puhti( user, password ):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect('puhti.csc.fi', username=user, password=password)
    sftp = client.open_sftp()
    return sftp

### Check file existence ###
def file_already_exists( filename ):
    return os.path.exists("./data/tribolium_eqtl/rds/"+filename)

##################
#      MAIN      #
##################

if __name__ == '__main__':
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments           #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(config["repository_path"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Load the list of phenotypes            #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load the list of phenotypes")
    phenotypes = load_phenotypes(config["dataset"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Connect to Puhti                       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Connect to Puhti")
    puhti = connect_to_puhti(config["user"], config["password"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 5) Download significant eQTL associations #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Download significant eQTL associations")
    rds_filename = config["dataset"]+"_significant.rds"
    csv_filename = config["dataset"]+"_significant.csv"
    puhti.get("/scratch/project_XXXXXXX/Tribolium_castaneum_eQTL/"+rds_filename, "./data/tribolium_eqtl/significant/"+rds_filename)
    puhti.get("/scratch/project_XXXXXXX/Tribolium_castaneum_eQTL/"+csv_filename, "./data/tribolium_eqtl/significant/"+csv_filename)

    puhti.close()

    print(">> Done.")

