#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 2_FilterGenotypes.py
# --------------------
# Filter SNPs based on genotype coverage.
# (LOCAL SCRIPT)
#
# 1) Filter genotypes based on F_MISSING,
# 2) Remove non-variant SNPs,
# 3) Remove LOWQUAL SNPs.
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
    parser.add_argument("--gatk-path", "-gatk-path", help="GATK path")
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--imputed", "-imputed", action="store_true")
    parser.add_argument("--suffix", "-suffix", help="Final suffix (also used to load filters)")
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

### Load filters ###
def load_filters( population, version, suffix ):
    ### Load filters ###
    filters  = {}
    filename = "./data/tribolium_filters/filters_"+population+"_"+version+"_"+suffix+".csv"
    f        = open(filename, "r")
    l        = f.readline()
    while l:
        l = l.strip("\n").split(";")
        filters[l[0]] = l[1]
        l = f.readline()
    f.close()
    ### Print filters for information ###
    print(" ---------------------------------------")
    print("| The following filters will be applied:")
    print(" ---------------------------------------")
    for item in filters.items():
        print("| • Filter '"+item[0]+"': "+item[1])
    print(" ---------------------------------------")
    return filters

### Mark as missing low read depth genotypes ###
def filter_genotype_DP( population, version, DP, in_suffix, out_suffix, delete_input ):
    input_vcf  = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    cmdline    = []
    cmdline.append("bcftools filter -i \"FORMAT/DP>="+str(DP)+"\"")
    cmdline.append("--set-GTs .")
    cmdline.append("-Ou "+input_vcf+" |")
    cmdline.append("bcftools view > "+output_vcf)
    os.system(" ".join(cmdline))
    assert_creation(output_vcf)
    if delete_input:
        delete_file(input_vcf)

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

### Filter based on genotype call rate (F_MISSING) ###
def filter_genotype_CR( population, version, F_MISSING, in_suffix, out_suffix, delete_input ):
    input_vcf  = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    cmdline    = []
    cmdline.append("bcftools view -i 'F_MISSING<"+str(F_MISSING)+"' "+input_vcf)
    cmdline.append("> "+output_vcf)
    os.system(" ".join(cmdline))
    assert_creation(output_vcf)
    if delete_input:
        delete_file(input_vcf)

### Run GATK VariantFiltration ###
def run_gatk_VariantFiltration( gatk_path, population, version, filters, in_suffix, out_suffix, delete_input ):
    input_vcf  = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    cmdline    = []
    cmdline.append(gatk_path+"gatk VariantFiltration")
    cmdline.append("-R ./data/tribolium_genome/Tribolium_castaneum_"+version+"/Tribolium_castaneum_"+version+".fna")
    cmdline.append("-V "+input_vcf)
    cmdline.append("-O "+output_vcf)
    for item in filters.items():
        if item[0] not in ["DP", "F_MISSING"]:
            cmdline.append("--filter-name "+item[0])
            cmdline.append("--filter-expression \""+item[1]+"\"")
    os.system(" ".join(cmdline))
    assert_creation(output_vcf)
    if delete_input:
        delete_file(input_vcf)

### Remove LOWQUAL SNPs ###
def remove_lowqual( population, version, in_suffix, out_suffix, delete_input ):
    input_vcf  = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    cmdline    = []
    cmdline.append("bcftools view -f PASS "+input_vcf+" > "+output_vcf)
    os.system(" ".join(cmdline))
    assert_creation(output_vcf)
    if delete_input:
        delete_file(input_vcf)

### Edit the VCF file ###
def edit_VCF_file( population, version, in_suffix, out_suffix, delete_input ):
    input_vcf  = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    f          = open(input_vcf, "r", encoding="ISO-8859-1")
    g          = open(output_vcf, "w", encoding="ISO-8859-1")
    l          = f.readline()
    while l:
        if "é" in l:
            l = l.replace("é", "e")
        if "à" in l:
            l = l.replace("à", "a")
        g.write(l)
        l = f.readline()
    f.close()
    g.close()
    if delete_input:
        delete_file(input_vcf)

### Compress a VCF file ###
def compress_VCF_file( population, version, suffix ):
    vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+suffix+".vcf"
    os.system("gzip -c "+vcf+" > "+vcf+".gz")

### Get the next suffix and the delete instruction ###
def next_suffix( current_suffix ):
    if current_suffix in ["raw_SNP", "imputed"]:
        return "1", False
    else:
        return str(int(current_suffix)+1), True


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
    SUFFIX_COUNTER = 1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Load filters                           #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load filters")
    filters = load_filters(config["population"], config["version"], config["suffix"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Define the initial suffix              #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    CURRENT_SUFFIX = "raw_SNP"
    DELETE         = False
    if config["imputed"]:
        CURRENT_SUFFIX = "imputed"

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Mark low coverage genotypes as missing #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if "DP" in filters.keys():
        print(">> Mark genotypes  with DP < "+filters["DP"]+" as missing")
        NEXT_SUFFIX, DELETE = next_suffix(CURRENT_SUFFIX)
        filter_genotype_DP(config["population"], config["version"], filters["DP"], CURRENT_SUFFIX, NEXT_SUFFIX, DELETE)
        CURRENT_SUFFIX = NEXT_SUFFIX
        print(">> Recalculate SNP tags")
        NEXT_SUFFIX, DELETE = next_suffix(CURRENT_SUFFIX)
        recalculate_SNP_tags(config["population"], config["version"], CURRENT_SUFFIX, NEXT_SUFFIX, DELETE)
        CURRENT_SUFFIX = NEXT_SUFFIX

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 5) Filter genotypes by call-rate          #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if "F_MISSING" in filters.keys():
        print(">> Filter genotypes with F_MISSING < "+str(filters["F_MISSING"]))
        NEXT_SUFFIX, DELETE = next_suffix(CURRENT_SUFFIX)
        filter_genotype_CR(config["population"], config["version"], filters["F_MISSING"], CURRENT_SUFFIX, NEXT_SUFFIX, DELETE)
        CURRENT_SUFFIX = NEXT_SUFFIX

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 6) Apply hard-filters                      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Apply hard-filters")
    NEXT_SUFFIX, DELETE = next_suffix(CURRENT_SUFFIX)
    run_gatk_VariantFiltration(config["gatk_path"], config["population"], config["version"], filters, CURRENT_SUFFIX, NEXT_SUFFIX, DELETE)
    CURRENT_SUFFIX = NEXT_SUFFIX

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 7) Remove low quality SNPs                 #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Remove low quality SNPs")
    NEXT_SUFFIX, DELETE = next_suffix(CURRENT_SUFFIX)
    remove_lowqual(config["population"], config["version"], CURRENT_SUFFIX, NEXT_SUFFIX, DELETE)
    CURRENT_SUFFIX = NEXT_SUFFIX

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 8) Edit VCF file                           #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Edit VCF file")
    edit_VCF_file(config["population"], config["version"], CURRENT_SUFFIX, config["suffix"], True)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 9) Compress resulting VCF file             #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Compress resulting VCF file")
    compress_VCF_file(config["population"], config["version"], config["suffix"])

    print(">> Done.")

