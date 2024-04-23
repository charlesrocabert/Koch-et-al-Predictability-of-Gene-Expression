#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 2_MergeFeatureCounts.py
# -----------------------
# Merge feature counts from every individual samples.
# (HPC SCRIPT --> run wrapper)
#
# 1) Import the list of samples and all read counts,
# 2) Merge read counts in a single file,
# 3) Export the resulting file to Allas.
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
    parser.add_argument("--suffix", "-suffix", help="Sample list filename suffix")
    parser.add_argument("--version", "-version", help="Reference genome version")
    args = parser.parse_args()
    return(vars(args))

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Import the list of samples ###
def import_sample_list( bucket, suffix, version ):
    os.system("a-get "+bucket+"/samples/samples_"+suffix+"_"+version+".csv")
    assert_creation("samples_"+suffix+"_"+version+".csv")

### Load the list of samples ###
def load_sample_list( suffix, version ):
    samples   = []
    file      = open("samples_"+suffix+"_"+version+".csv", "r")
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

### Import all the individual read counts ###
def import_read_counts_folder( bucket ):
    os.system("mkdir read_counts")
    os.system("rclone copyto allas:"+bucket+"/read_counts read_counts/.")

### Merge read counts ###
def merge_read_counts( samples, bucket, suffix, version ):
    #---------------------------------------#
    # 1) Parse individual read counts       #
    #---------------------------------------#
    read_counts = {}
    summary     = {}
    first_pass  = True
    for sample in samples[1:]:
        sample_name = sample["sample"]
        print("   > Parsing sample "+sample_name)
        f = open("read_counts/"+sample_name+".txt", "r")
        l = f.readline()
        l = f.readline()
        l = f.readline()
        while l:
            l       = l.strip("\n").split("\t")
            gene_id = l[0]
            count   = int(l[6])
            if first_pass:
                assert gene_id not in read_counts.keys()
                read_counts[gene_id] = {}
                read_counts[gene_id][sample_name] = count
            else:
                assert gene_id in read_counts.keys()
                assert sample_name not in read_counts[gene_id].keys()
                read_counts[gene_id][sample_name] = count
            l  = f.readline()
        f.close()
        first_pass = False
    #---------------------------------------#
    # 2) Merge read counts in a single file #
    #---------------------------------------#
    g      = open("Tribolium_castaneum_"+suffix+"_"+version+"_read_counts.txt", "w")
    header = "gene_id"
    for sample in samples[1:]:
        header += "\t"+sample["sample"]
    g.write(header+"\n")
    for gene_id in read_counts.keys():
        l = gene_id
        for sample in samples[1:]:
            sample_name = sample["sample"]
            l += "\t"+str(read_counts[gene_id][sample_name])
        g.write(l+"\n")
    g.close()

### Export merged read counts file ###
def export_merged_read_counts( bucket, suffix, version ):
    txt_filename = "Tribolium_castaneum_"+suffix+"_"+version+"_read_counts.txt"
    os.system("rclone copyto "+txt_filename+" allas:"+bucket+"/read_counts/"+txt_filename)


##################
#      MAIN      #
##################

if __name__ == '__main__':
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(os.environ['DATADIR'])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Import the list of samples        #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Import the list of samples")
    import_sample_list(config["bucket"], config["suffix"], config["version"])
    samples = load_sample_list(config["suffix"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Import all individual read counts #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Import all individual read counts")
    import_read_counts_folder(config["bucket"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Merge and export read counts      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    ### 4.1) Merge the read counts ###
    print(">> Merge the read counts")
    merge_read_counts(samples, config["bucket"], config["suffix"], config["version"])

    ### 4.2) Export merged read counts file ###
    print(">> Export merged read counts file")
    export_merged_read_counts(config["bucket"], config["suffix"], config["version"])

    print(">> Done.")

