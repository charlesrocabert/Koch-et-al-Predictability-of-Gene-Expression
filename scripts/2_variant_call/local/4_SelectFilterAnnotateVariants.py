#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 4_SelectFilterAnnotateVariants.py
# ---------------------------------
# Select bi-allelic SNP variants from the original VCF file.
# (LOCAL SCRIPT)
#
# 1) Select bi-allelic SNP variants
# 2) Tag DP=0 genotypes as missing ('./.'),
# 3) Filter out low quality SNPs,
# 4) Annotate SNPs,
# 5) Add SNP unique identifiers.
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
    parser.add_argument("--gatk-path", "-gatk-path", help="GATK folder path")
    parser.add_argument("--snpeff-path", "-snpeff-path", help="snpeff .jar file path")
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
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

### Run GATK SelectVariants ###
def run_gatk_SelectVariants( gatk_path, population, version, out_suffix ):
    input_vcf  = "./data/tribolium_vcf/Tribolium_castaneum_"+population+"_"+version+".vcf"
    output_vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    cmdline    = []
    cmdline.append(gatk_path+"/gatk SelectVariants")
    cmdline.append("-R ./data/tribolium_genome/Tribolium_castaneum_"+version+"/Tribolium_castaneum_"+version+".fna")
    cmdline.append("-V "+input_vcf)
    cmdline.append("-O "+output_vcf)
    cmdline.append("--select-type-to-include SNP")
    cmdline.append("--restrict-alleles-to BIALLELIC")
    cmdline.append("--exclude-non-variants")
    os.system(" ".join(cmdline))
    assert_creation(output_vcf)

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

### Annotate the VCF file with snpEff ###
def annotate( snpeff_path, population, version, in_suffix, out_suffix, delete_input ):
    #------------------------------#
    # 1) Build the snpEff database #
    #------------------------------#
    os.chdir("./data/tribolium_snp")
    os.system("mkdir ./data")
    os.system("mkdir ./data/tribolium")
    f = open("snpEff.config", "w")
    f.write("tribolium.genome : tribolium\n")
    f.write("data.dir = ./data/\n")
    f.close()
    os.system("cp ../tribolium_genome/Tribolium_castaneum_"+version+"/Tribolium_castaneum_"+version+".gff3 data/tribolium/genes.gff")
    os.system("cp ../tribolium_genome/Tribolium_castaneum_"+version+"/Tribolium_castaneum_"+version+".fna data/tribolium/sequences.fa")
    os.system("cp "+snpeff_path+" snpEff.jar")
    os.system("java -jar snpEff.jar build -c snpEff.config -gff3 -v tribolium")
    #------------------------------#
    # 2) Run the annotation        #
    #------------------------------#
    input_vcf  = "Tribolium_castaneum_"+population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    cmdline    = []
    cmdline.append("java -jar snpEff.jar -nodownload -c snpEff.config tribolium "+input_vcf+" > "+output_vcf)
    os.system(" ".join(cmdline))
    assert_creation(output_vcf)
    if delete_input:
        delete_file(input_vcf)
    delete_file("snpEff.config")
    delete_file("snpEff.jar")
    os.chdir("../..")

### Add variant IDs to the VCF file ###
def add_variant_identifiers( population, version, in_suffix, out_suffix, delete_input ):
    input_vcf  = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+in_suffix+".vcf"
    output_vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+out_suffix+".vcf"
    f          = open(input_vcf, "r", encoding="ISO-8859-1")
    g          = open(output_vcf, "w", encoding="ISO-8859-1")
    ### Skip the header ###
    l = f.readline()
    while "#CHROM" not in l:
        g.write(l)
        l = f.readline()
    g.write(l)
    ### Parse the file ###
    counter = 1
    l       = f.readline()
    while l:
        if counter%10000 == 0:
            print("   > "+str(counter)+" variants parsed")
        l = l.strip("\n").split("\t")
        l[2] = l[0]+"-"+l[1]
        g.write("\t".join(l)+"\n")
        l = f.readline()
        counter += 1
    f.close()
    g.close()
    print("   > "+str(counter)+" SNPs have been parsed.")
    if delete_input:
        delete_file(input_vcf)

### Compress a VCF file ###
def compress_VCF_file( population, version, suffix ):
    vcf = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+suffix+".vcf"
    os.system("gzip -c "+vcf+" > "+vcf+".gz")


##################
#      MAIN      #
##################

if __name__ == '__main__':

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config         = parse_arguments()
    SUFFIX_COUNTER = 1
    os.chdir(config["repository_path"])
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Load filters                   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load filters")
    filters = load_filters(config["population"], config["version"], config["suffix"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Run GATK SelectVariants        #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Select bi-allelic variants")
    run_gatk_SelectVariants(config["gatk_path"], config["population"], config["version"], str(SUFFIX_COUNTER))
    SUFFIX_COUNTER += 1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Mark DP=0 genotypes as missing #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Mark DP=0 genotypes as missing")
    filter_genotype_DP(config["population"], config["version"], filters["DP"], str(SUFFIX_COUNTER-1), str(SUFFIX_COUNTER), True)
    SUFFIX_COUNTER += 1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 5) Recalculate SNP tags           #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Recalculate SNP tags")
    recalculate_SNP_tags(config["population"], config["version"], str(SUFFIX_COUNTER-1), str(SUFFIX_COUNTER), True)
    SUFFIX_COUNTER += 1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 6) Apply hard-filters             #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Apply hard-filters")
    run_gatk_VariantFiltration(config["gatk_path"], config["population"], config["version"], filters, str(SUFFIX_COUNTER-1), str(SUFFIX_COUNTER), True)
    SUFFIX_COUNTER += 1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 7) Remove low quality SNPs        #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Remove low quality SNPs")
    remove_lowqual(config["population"], config["version"], str(SUFFIX_COUNTER-1), str(SUFFIX_COUNTER), True)
    SUFFIX_COUNTER += 1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 8) Edit VCF file                  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Edit VCF file")
    edit_VCF_file(config["population"], config["version"], str(SUFFIX_COUNTER-1), str(SUFFIX_COUNTER), True)
    SUFFIX_COUNTER += 1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 9) Annotate SNPs                  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Annotate SNPs")
    annotate(config["snpeff_path"], config["population"], config["version"], str(SUFFIX_COUNTER-1), str(SUFFIX_COUNTER), True)
    SUFFIX_COUNTER += 1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 10) Add SNP unique identifiers    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Add SNP unique identifiers")
    add_variant_identifiers(config["population"], config["version"], str(SUFFIX_COUNTER-1), config["suffix"], True)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 11) Compress resulting VCF file   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Compress resulting VCF file")
    compress_VCF_file(config["population"], config["version"], config["suffix"])

    print(">> Done.")

