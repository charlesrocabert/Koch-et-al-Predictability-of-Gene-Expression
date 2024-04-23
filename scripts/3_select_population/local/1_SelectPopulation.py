#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_SelectPopulation.py
# ---------------------
# Select a population of samples based on a samples list.
# (LOCAL SCRIPT)
#
# 1) Select the subset of samples in the VCF file,
# OR
# 2) Select the subset of samples in the read counts file.
#***************************************************************************

import os
import sys
import csv
import time
import argparse

### Parse command line arguments ###
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-path", "-repository-path", help="Repository path")
    parser.add_argument("--in-population", "-in-population", help="Input population")
    parser.add_argument("--in-suffix", "-in-suffix", help="Input population suffix")
    parser.add_argument("--out-population", "-out-population", help="Output population")
    parser.add_argument("--out-suffix", "-out-suffix", help="Output population suffix")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--task", "-task", choices=['VCF', 'COUNTS', 'BOTH'], help="Which task(s) to run")
    args   = parser.parse_args()
    config = vars(args)
    assert config["in_population"]+"_"+config["in_suffix"] != config["out_population"]+"_"+config["out_suffix"], "Input and output must be different. Exit."
    return(config)

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Delete a file ###
def delete_file( filename ):
    os.system("rm "+filename)
    assert_deletion(filename)

### Load the list of samples ###
def load_samples_list( out_population, version ):
    samples   = []
    file      = open("./data/tribolium_bam/samples_"+out_population+"_"+version+".csv", "r")
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

### Select the sub-population in the VCF file ###
def select_sub_population_in_VCF( samples, in_population, in_suffix, out_population, out_suffix, version ):
    ### Generate samples list ###
    global_params = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    sample_list   = []
    for sample in samples[1:]:
        sample_list.append(sample["sample"])
    ### Open input and output VCF files ###
    input_vcf  = "./data/tribolium_snp/Tribolium_castaneum_"+in_population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "./data/tribolium_snp/Tribolium_castaneum_"+out_population+"_"+version+"_"+out_suffix+".vcf"
    f          = open(input_vcf, "r", encoding="ISO-8859-1")
    g          = open(output_vcf, "w", encoding="ISO-8859-1")
    ### Skip the header ###
    l = f.readline()
    while "#CHROM" not in l:
        g.write(l)
        l = f.readline()
    ### Parse the header ###
    header          = l.strip("\n").split("\t")
    samples_to_keep = []
    for i in range(len(header)):
        if header[i] in global_params or header[i] in sample_list:
            samples_to_keep.append(i)
    counter = 1
    while l:
        if counter%10000 == 0:
            print("   > "+str(counter)+" variants parsed")
        l = l.strip("\n").split("\t")
        if l[0] != "Unknown":
            line = ""
            for i in range(len(l)):
                if i in samples_to_keep:
                    line += l[i]+"\t"
            g.write(line.strip("\t")+"\n")
        l = f.readline()
        counter += 1
    f.close()
    g.close()
    os.system("gzip -c "+output_vcf+" > "+output_vcf+".gz")

### Recalculate SNP tags ###
def recalculate_SNP_tags( population, version, in_suffix, out_suffix, delete_input ):
    input_vcf  = "data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    cmdline    = []
    cmdline.append("bcftools +fill-tags "+input_vcf+" |")
    cmdline.append("bcftools view > "+output_vcf)
    os.system(" ".join(cmdline))
    assert_creation(output_vcf)
    if delete_input:
        delete_file(input_vcf)

### Rename the VCF file ###
def rename_vcf_file( population, version, in_suffix, out_suffix ):
    input_vcf  = "data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    os.system("mv "+input_vcf+" "+output_vcf)
    assert_deletion(input_vcf)
    assert_creation(output_vcf)

### Select the sub-population in the read counts file ###
def select_sub_population_in_read_counts( samples, in_pop, in_suffix, out_pop, out_suffix, version ):
    ### Generate samples list ###
    global_params = ["gene_id"]
    sample_list   = []
    for sample in samples[1:]:
        sample_list.append(sample["sample"])
    ### Open input and output VCF files ###
    input_txt  = "./data/tribolium_counts/Tribolium_castaneum_"+in_pop+"_"+version+"_"+in_suffix+".txt"
    output_txt = "./data/tribolium_counts/Tribolium_castaneum_"+out_pop+"_"+version+"_"+out_suffix+".txt"
    f = open(input_txt, "r")
    g = open(output_txt, "w")
    ### Parse the header ###
    l               = f.readline()
    header          = l.strip("\n").split("\t")
    samples_to_keep = []
    for i in range(len(header)):
        if header[i] in global_params or header[i] in sample_list:
            samples_to_keep.append(i)
    counter = 1
    while l:
        if counter%1000 == 0:
            print("   > "+str(counter)+" transcripts parsed")
        l    = l.strip("\n").split("\t")
        line = ""
        for i in range(len(l)):
            if i in samples_to_keep:
                line += l[i]+"\t"
        g.write(line.strip("\t")+"\n")
        l = f.readline()
        counter += 1
    f.close()
    g.close()
    os.system("gzip -c "+output_txt+" > "+output_txt+".gz")


##################
#      MAIN      #
##################

if __name__ == '__main__':
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments                  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(config["repository_path"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Import the list of samples and                #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Import the list of samples")
    samples  = load_samples_list(config["out_population"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Select the sub-population in VCF file         #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if config["task"] in ["VCF", "BOTH"]:
        print(">> Select the sub-population in VCF file")
        select_sub_population_in_VCF(samples, config["in_population"], config["in_suffix"], config["out_population"], config["out_suffix"], config["version"])
        recalculate_SNP_tags(config["out_population"], config["version"], config["out_suffix"], config["out_suffix"]+"_retag", True)
        rename_vcf_file(config["out_population"], config["version"], config["out_suffix"]+"_retag", config["out_suffix"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Select the sub-population in read counts file #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if config["task"] in ["COUNTS", "BOTH"]:
        print(">> Select the sub-population in read counts file")
        select_sub_population_in_read_counts(samples, config["in_population"], config["in_suffix"], config["out_population"], config["out_suffix"], config["version"])

    print(">> Done.")

