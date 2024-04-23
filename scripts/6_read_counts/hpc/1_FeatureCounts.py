#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_FeatureCounts.py
# ------------------
# Calculate feature counts for every individual samples.
# (HPC SCRIPT --> run wrapper)
#
# 1) Import subread package and compile it,
# 2) Import reference genome annotation,
# 3) Import the list of samples,
# 4) For each BAM file, if the read counts file does not exist:
#   4.1) Run subread FeatureCounts,
#   4.2) Export resulting files to Allas.
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

### Import subread package and compile it ###
def import_subread_package( bucket ):
    os.system("a-get "+bucket+"/read_counts/subread_package/subread-2.0.3-source.tar.gz")
    os.system("tar -xf subread-2.0.3-source.tar.gz")
    os.chdir("subread-2.0.3-source/src/")
    os.system("make -f Makefile.Linux")
    os.chdir("../..")
    assert_creation("subread-2.0.3-source/bin/featureCounts")

### Import reference genome annotation ###
def import_reference_genome_annotation( bucket, version ):
    genome_annotation = "Tribolium_castaneum_"+version+".gtf"
    os.system("a-get "+bucket+"/reference_genome/"+genome_annotation+".zip")
    assert_creation(genome_annotation+".zip")
    os.system("unzip "+genome_annotation+".zip")
    assert_creation(genome_annotation)

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

### Import BAM file ###
def import_bam( input_bam, bucket ):
    os.system("a-get "+bucket+"/bam/"+input_bam)
    assert_creation(input_bam)

### Check if the read counts have been calculated ###
def readcounts_file_exists( bucket, sample_name ):
    cmdline = "a-list "+bucket+"/read_counts/"+sample_name+".txt"
    process = subprocess.Popen([cmdline], stdout=subprocess.PIPE, shell=True)
    output  = process.stdout.read().decode('utf8').strip("\n")
    if len(output) == 0:
        return False
    else:
        return True

### Run subread featureCounts ###
def run_subread_featureCounts( sample_name, input_bam, version ):
    genome_annotation = "Tribolium_castaneum_"+version+".gtf"
    cmdline           = []
    cmdline.append("subread-2.0.3-source/bin/featureCounts")
    #############
    cmdline.append("-s 2")
    cmdline.append("--minOverlap 10")
    cmdline.append("-Q 10")
    cmdline.append("-O")
    cmdline.append("-M")
    cmdline.append("--primary")
    #############
    cmdline.append("-a Tribolium_castaneum_"+version+".gtf")
    cmdline.append("-o "+sample_name+".txt")
    cmdline.append(" "+input_bam)
    os.system(" ".join(cmdline))
    assert_creation(sample_name+".txt")
    assert_creation(sample_name+".txt.summary")

### Export read counts files ###
def export_read_counts( sample_name, bucket ):
    os.system("rclone copyto "+sample_name+".txt allas:"+bucket+"/read_counts/"+sample_name+".txt")
    os.system("rclone copyto "+sample_name+".txt.summary allas:"+bucket+"/read_counts/"+sample_name+".txt.summary")


##################
#      MAIN      #
##################

if __name__ == '__main__':
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments                  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(os.environ['DATADIR'])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Import needed package and datasets            #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    ### 2.1) Import subread package and compile it ###
    print(">> Import subread package and compile it")
    import_subread_package(config["bucket"])

    ### 2.2) Import reference genome annotation ###
    print(">> Import reference genome annotation")
    import_reference_genome_annotation(config["bucket"], config["version"])

    ### 2.3) Import the list of samples ###
    print(">> Import the list of samples")
    import_sample_list(config["bucket"], config["suffix"], config["version"])
    samples = load_sample_list(config["suffix"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Compute subread featureCounts for each sample #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    for sample in samples[1:]:
        sample_name = sample["sample"]
        if not readcounts_file_exists(config["bucket"], sample_name):
            print(">> Read counts file does not exists, run the pipeline")

            ### 3.1) Import BAM file ###
            print("   > Import the BAM file")
            import_bam(sample_name+"_FC.bam", config["bucket"])

            ### 3.2) Run subread featureCounts ###
            print("   > Run subread featureCounts")
            run_subread_featureCounts(sample_name, sample_name+"_FC.bam", config["version"])

            ### 3.3) Export read counts files ###
            print("   > Export read counts file")
            export_read_counts(sample_name, config["bucket"])

    print(">> Done.")

