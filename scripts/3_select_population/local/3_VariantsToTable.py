#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 3_VariantsToTable.py
# --------------------
# Extract VCF statistics into a text table.
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
    parser.add_argument("--gatk-path", "-gatk-path", help="GATK path")
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--suffix", "-suffix", help="Final suffix")
    parser.add_argument("--annotation", "-annotation", action="store_true")
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

### Run GATK VariantsToTable ###
def run_gatk_VariantsToTable( gatk_path, population, version, suffix ):
    input   = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+suffix+".vcf"
    output  = "./data/tribolium_snp/snp_table_"+population+"_"+version+"_"+suffix+".csv"
    stats   = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "AC", "AF", "AN", "BaseQRankSum", "DP", "ExcessHet", "FS", "InbreedingCoeff", "MLEAC", "MLEAF", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR", "HET", "HOM-REF", "HOM-VAR", "NO-CALL", "VAR", "NSAMPLES", "NCALLED", "HWE", "MAF", "ANN"]
    cmdline = []
    cmdline.append(gatk_path+"/gatk VariantsToTable")
    cmdline.append("-V "+input)
    cmdline.append("-O "+output)
    for elmt in stats:
        cmdline.append("-F "+elmt)
    os.system(" ".join(cmdline))
    assert_creation(output)

### Correct GATK header bug ###
def correct_table( population, version, suffix ):
    tablename = "./data/tribolium_snp/snp_table_"+population+"_"+version+"_"+suffix+".csv"
    f = open(tablename, "r")
    g = open(tablename+"_edit", "w")
    l = f.readline()
    g.write(l.strip("\t\n")+"\n")
    l = f.readline()
    while l:
        g.write(l)
        l = f.readline()
    f.close()
    g.close()
    os.system("mv "+tablename+"_edit"+" "+tablename)

### Parse annotations in the VCF table ###
def parse_annotations( population, version, suffix ):
    tablename  = "./data/tribolium_snp/snp_table_"+population+"_"+version+"_"+suffix+".csv"
    ann_list   = ["Annotation", "Putative_impact", "Gene_name", "Gene_id", "Feature_type", "Feature_id", "Transcript_biotype", "Rank_total", "HGVS_c", "HGVS_p", "cDNA_position", "CDS_position", "Protein_position", "Distance_to_feature", "Errors"]
    f          = open(tablename, "r")
    g          = open(tablename+"_edit", "w")
    header     = f.readline()
    new_header = header.strip("\n").split("\t")
    new_header = "\t".join(new_header[0:(len(new_header)-1)])+"\t"+"\t".join(ann_list)+"\n"
    ann_index  = len(header.strip("\n").split("\t"))-1
    g.write(new_header)
    l = f.readline()
    while l:
        new_line = l.strip("\n").split("\t")
        new_line = "\t".join(new_line[0:(len(new_line)-1)])
        ann      = l.strip("\n").split("\t")[ann_index].split(",")[0].split("|")
        for i in range(len(ann_list)):
            if ann_list[i] == "Feature_id":# and "transcript:" in ann[(i+1)]:
                #new_line += "\t"+ann[(i+1)].split("-")[0].split(":")[1]
                new_line += "\t"+ann[(i+1)].strip("transcript:").strip("-RA")
            else:
                if ann[(i+1)] == "":
                    new_line += "\t-"
                else:
                    new_line += "\t"+ann[(i+1)]
        g.write(new_line+"\n")
        l = f.readline()
    f.close()
    g.close()
    os.system("mv "+tablename+"_edit "+tablename)


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
    # 2) Run GATK SelectVariants      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Extract VCF statistics into a text table")
    run_gatk_VariantsToTable(config["gatk_path"], config["population"], config["version"], config["suffix"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Remove GATK header bug       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Remove GATK header bug")
    correct_table(config["population"], config["version"], config["suffix"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Parse annotations            #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if config["annotation"]:
        print(">> Parse snpEff annotations")
        parse_annotations(config["population"], config["version"], config["suffix"])

    print(">> Done.")

