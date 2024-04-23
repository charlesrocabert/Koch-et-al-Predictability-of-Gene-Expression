#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 2_UploadGemmaFiles.py
# ---------------------
# Upload gemma input files to a distant server.
# (LOCAL SCRIPT)
#***************************************************************************

import os
import sys
import csv
import time
import argparse
import subprocess
import swiftclient

### Parse command line arguments ###
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-path", "-repository-path", help="Path to the repository")
    parser.add_argument("--user", "-user", help="Allas user name")
    parser.add_argument("--password", "-password", help="Allas password")
    parser.add_argument("--project", "-project", help="Allas project")
    parser.add_argument("--bucket", "-bucket", help="Allas bucket name")
    args = parser.parse_args()
    return(vars(args))

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Connect to Allas ###
def connect_to_allas( user, password, project ):
    _authurl      = 'https://pouta.csc.fi:5001/v3'
    _auth_version = '3'
    _user         = user
    _key          = password
    _os_options   = {
        'user_domain_name': 'Default',
        'project_domain_name': 'Default',
        'project_name': project
    }
    conn = swiftclient.Connection(
        authurl=_authurl,
        user=_user,
        key=_key,
        os_options=_os_options,
        auth_version=_auth_version
    )
    return conn

### Put the file on Allas ###
def put_on_allas( allas_connection, bucket, filename ):
    open_key = "r"
    if filename.endswith(".bed") or filename.endswith(".bim"):
        open_key = "rb"
    with open("./data/tribolium_eqtl/gemma/"+filename, open_key) as f:
        allas_connection.put_object(bucket, "eqtl/"+filename, contents=f.read())


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
    # 2) Connect to Allas             #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Connect to Allas")
    allas = connect_to_allas(config["user"], config["password"], config["project"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Export input files to Allas  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Export input files to Allas")
    ENVIRONMENT = ["CT", "HD"]
    IMPUTATION  = ["_imputed"]
    PHENOTYPE   = ["EXPRESSION", "PLASTICITY", "FITNESS", "NOISE"]
    EXTENSION   = ["bed", "bim", "fam", "pheno"]
    for env in ENVIRONMENT:
        for imp in IMPUTATION:
            for phe in PHENOTYPE:
                for ext in EXTENSION:
                    if not (env=="CT" and phe in ["PLASTICITY","NOISE"]):
                        filename = env+"_G1_Tcas3.30"+imp+"_"+phe+"."+ext
                        print("   > Dealing with file "+filename)
                        put_on_allas(allas, config["bucket"], filename)

    print(">> Done.")

