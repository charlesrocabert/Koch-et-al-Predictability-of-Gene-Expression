#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 2_HaplotypeCaller.py")
# --------------------")
# Run the complete pipeline for per-sample variant call.
# (HPC SCRIPT --> array wrapper)
#
# 1) Haplotype-caller (on the marked duplicates BAM file):
#   1.1) Import the BAM file from Allas
#   1.2) Import the reference genome and generate indexes,
#   1.3) Generate BAI index file,
#   1.4) Run GATK HaplotypeCaller. This task can take several hours,
#   1.5) Export GVCF files to Allas.
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
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--sindex", "-sindex", type=int, default=0, help="Index of the sample in the list of samples")
    args = parser.parse_args()
    return(vars(args))

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Import the list of samples ###
def import_sample_list( bucket, population, version ):
    os.system("a-get "+bucket+"/samples/samples_"+population+"_"+version+".csv")
    assert_creation("samples_"+population+"_"+version+".csv")

### Load the list of samples ###
def load_sample_list( population, version ):
    samples   = []
    file      = open("samples_"+population+"_"+version+".csv", "r")
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
def import_bam( input_bam, output_bam, bucket ):
    os.system("a-get "+bucket+"/bam/"+input_bam)
    assert_creation(input_bam)
    if input_bam != output_bam:
        os.system("mv "+input_bam+" "+output_bam)
        assert_deletion(input_bam)
        assert_creation(output_bam)

### Import reference genome ###
def import_reference_genome( bucket, version ):
    reference_genome = "Tribolium_castaneum_"+version+".fna"
    os.system("a-get "+bucket+"/reference_genome/"+reference_genome+".zip")
    assert_creation(reference_genome+".zip")
    os.system("unzip "+reference_genome+".zip")
    assert_creation(reference_genome)

### Create reference genome indices ###
def create_reference_genome_indices( version ):
    reference_genome = "Tribolium_castaneum_"+version+".fna"
    os.system("samtools faidx "+reference_genome)
    os.system("gatk CreateSequenceDictionary -R "+reference_genome)
    assert_creation(reference_genome+".fai")
    assert_creation(reference_genome.strip(".fna")+".dict")

### Create BAM file index ###
def create_bam_index( input_bam ):
    os.system("samtools index "+input_bam)
    assert_creation(input_bam+".bai")

### Run GATK HaplotypeCaller ###
def run_GATK_HaplotypeCaller( input_bam, output_gvcf, version ):
    reference_genome = "Tribolium_castaneum_"+version+".fna"
    cmdline          = []
    cmdline.append("gatk HaplotypeCaller")
    cmdline.append("-R "+reference_genome)
    cmdline.append("-I "+input_bam)
    cmdline.append("-O "+output_gvcf)
    cmdline.append("-ERC GVCF")
    cmdline.append("--pcr-indel-model NONE")
    cmdline.append("--dont-use-soft-clipped-bases")
    os.system(" ".join(cmdline))
    assert_creation(output_gvcf)
    assert_creation(output_gvcf+".tbi")

### Export BAM file to Allas ###
def export_per_sample_gvcf( input_gvcf, bucket ):
    os.system("rclone copyto "+input_gvcf+" allas:"+bucket+"/gvcf/"+input_gvcf)
    os.system("rclone copyto "+input_gvcf+".tbi allas:"+bucket+"/gvcf/"+input_gvcf+".tbi")


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
    # 2) Import the list of samples and get the sample #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Import the list of samples")
    import_sample_list(config["bucket"], config["population"], config["version"])
    samples     = load_sample_list(config["population"], config["version"])
    sample      = samples[config["sindex"]]
    sample_name = sample["sample"]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Run GATK HaplotypeCaller                      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    ### 3.1) Import data ###
    print(">> Import BAM file")
    import_bam(sample_name+".bam", sample_name+".bam", config["bucket"])

    print(">> Import reference genome and create indices")
    import_reference_genome(config["bucket"], config["version"])
    create_reference_genome_indices(config["version"])

    ### 3.2) Create BAM index ###
    print(">> Create BAM index")
    create_bam_index(sample_name+".bam")

    ### 3.3) Run GATK HaplotypeCaller in GVCF mode ###
    print(">> Run GATK HaplotypeCaller")
    run_GATK_HaplotypeCaller(sample_name+".bam", sample_name+".g.vcf.gz", config["version"])

    ### 3.4) Export the resulting GVCF file ###
    print(">> Export GVCF")
    export_per_sample_gvcf(sample_name+".g.vcf.gz", config["bucket"])

    print(">> Done.")

