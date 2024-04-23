#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 2_ImputeGenotypes.py
# --------------------
# Impute genotypes of a VCF file with Beagle.
# (LOCAL SCRIPT)
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
    parser.add_argument("--repository-path", "-repository-path", help="Repository path")
    parser.add_argument("--beagle", "-beagle", help="Beagle-5.4 path")
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--suffix", "-suffix", help="Suffix")
    parser.add_argument("--map", "-map", action="store_true")
    args = parser.parse_args()
    return(vars(args))

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

### Run Beagle for imputation ###
def run_Beagle_imputation( beagle_path, population, version, suffix, map ):
    vcf_file = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+suffix+".vcf"
    map_file = ""
    if map:
        map_file = "./data/tribolium_ld/Tribolium_castaneum_"+version+"_plink_map.txt"
    cmdline = []
    cmdline.append("java -jar "+beagle_path)
    cmdline.append("gt="+vcf_file)
    cmdline.append("out=./data/tribolium_snp/imputed_"+population+"_"+version+"_"+suffix)
    if map:
        cmdline.append("map="+map_file)
    os.system(" ".join(cmdline))
    assert_creation("./data/tribolium_snp/imputed_"+population+"_"+version+"_"+suffix+".vcf.gz")
    os.system("gunzip ./data/tribolium_snp/imputed_"+population+"_"+version+"_"+suffix+".vcf.gz")
    assert_deletion("./data/tribolium_snp/imputed_"+population+"_"+version+"_"+suffix+".vcf.gz")

### Load the genotypes from a VCF file ###
def load_genotypes( population, version, suffix ):
    vcf       = "./data/tribolium_snp/imputed_"+population+"_"+version+"_"+suffix+".vcf"
    genotypes = {"header":"", "missing":[]}
    f         = open(vcf, "r", encoding="ISO-8859-1")
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

### Write the merged VCF file ###
def merge_and_write_VCF_file( imputed_genotypes, population, version, suffix ):
    input_vcf     = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+suffix+".vcf"
    output_vcf    = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_merged.vcf"
    f             = open(input_vcf, "r", encoding="ISO-8859-1")
    g             = open(output_vcf, "w", encoding="ISO-8859-1")
    global_params = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    ### Save the header lines ###
    l = f.readline()
    while "#CHROM" not in l:
        g.write(l)
        l = f.readline()
    ### Parse and save the header ###
    g.write(l)
    header  = l.strip("\n").split("\t")
    l       = f.readline()
    counter = 1
    while l:
        if counter%10000 == 0:
            print("   > "+str(counter)+" variants parsed")
        l  = l.strip("\n").split("\t")
        ID = l[2]
        assert ID in imputed_genotypes.keys()
        line = ""
        for i in range(len(header)):
            if header[i] in global_params:
                line += l[i]+"\t"
            else:
                assert header[i] in imputed_genotypes[ID].keys()
                old_geno  = l[i].split(":")[0]
                new_geno  = imputed_genotypes[ID][header[i]]
                info      = l[i].strip(old_geno+":")
                line     += new_geno+":"+info+"\t"
        g.write(line.strip("\t")+"\n")
        l = f.readline()
        counter += 1
    f.close()
    g.close()

### Recalculate SNP tags ###
def recalculate_SNP_tags( population, version, in_suffix, out_suffix, delete_input ):
    input_vcf  = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    cmdline    = []
    cmdline.append("bcftools +fill-tags "+input_vcf+" |")
    cmdline.append("bcftools view > "+output_vcf)
    os.system(" ".join(cmdline))
    assert_creation(output_vcf)
    if delete_input:
        delete_file(input_vcf)

### Delete intermediate files ###
def delete_intermediate_files( population, version, suffix ):
    os.system("rm data/tribolium_snp/imputed_"+population+"_"+version+"_"+suffix+".vcf")
    os.system("rm data/tribolium_snp/imputed_"+population+"_"+version+"_"+suffix+".vcf.gz")
    os.system("rm data/tribolium_snp/imputed_"+population+"_"+version+"_"+suffix+".log")

### Compress a VCF file ###
def compress_VCF_file( population, version, suffix ):
    vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+suffix+".vcf"
    os.system("gzip -c "+vcf+" > "+vcf+".gz")


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
    # 2) Impute genotypes             #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Impute genotypes with Beagle")
    run_Beagle_imputation(config["beagle"], config["population"], config["version"], config["suffix"], config["map"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Load imputed genotypes       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load imputed genotypes")
    imputed_genotypes = load_genotypes(config["population"], config["version"], config["suffix"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Merge and write genotypes    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Merge and write genotypes in a single VCF file")
    merge_and_write_VCF_file(imputed_genotypes, config["population"], config["version"], config["suffix"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 5) Recalculate SNP tags         #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Recalculate SNP tags")
    recalculate_SNP_tags(config["population"], config["version"], "merged", "imputed", True)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 6) Delete intermediate files    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Delete intermediate files")
    delete_intermediate_files(config["population"], config["version"], config["suffix"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 7) Compress resulting VCF file  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Compress resulting VCF file")
    compress_VCF_file(config["population"], config["version"], "imputed")

    print(">> Done.")

