#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 3_JointCaller.py
# ----------------
# Run the complete pipeline for the join call.
# (HPC SCRIPT --> run wrapper)
#
# 1) Generate GenomicsDB database:
#   1.1) Import the list of samples,
#   1.2) Import all GVCFs from Allas,
#   1.3) Import the reference genome from Allas and compute index files,
#   1.4) Generate the sample map,
#   1.5) Generate the interval list,
#   1.6) Run GATK GenomicsDBImport,
#   1.7) Export the database to the scratch.
# 2) Joint-call cohort:
#   2.1) Import the consolidated database from the scratch if needed,
#   2.2) Run GATK GenotypeGVCFs,
#   2.3) Export the joint-call file (VCF) to the scratch.
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
    parser.add_argument("--task", "-task", choices=['GENOMICS_DB', 'GENOTYPE_GVCFS', 'BOTH'], help="Which task(s) to run")
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

### Generate the GVCF sample map ###
def generate_gvcf_sample_map( samples ):
    f = open("gvcf_map.csv", "w")
    for sample in samples[1:]:
        f.write(sample["sample"]+"\tgvcf/"+sample["sample"]+".g.vcf.gz\n")
    f.close()

### Import all the per-sample GVCF (and their indices) ###
def import_gvcf_folder( bucket ):
    os.system("mkdir gvcf")
    os.system("rclone copyto allas:"+bucket+"/gvcf gvcf/.")

### Check if all per-sample GVCF files have been imported ###
def check_gvcf_files( samples ):
    for sample in samples[1:]:
        assert_creation("gvcf/"+sample["sample"]+".g.vcf.gz")
        assert_creation("gvcf/"+sample["sample"]+".g.vcf.gz.tbi")

### Delete the per-sample GVCFs folder ###
def delete_gvcf_folder():
    os.system("rm -rf gvcf")

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

### Run GATK GenomicsDBImport ###
def run_gatk_GenomicsDBImport( version ):
    f = open("interval.list", "w")
    #f.write("ChLGX\nChLG2\nChLG3\nChLG4\nChLG5\nChLG6\nChLG7\nChLG8\nChLG9\nChLG10\nUnknown\n")
    f.write("ChLGX\nChLG2\nChLG3\nChLG4\nChLG5\nChLG6\nChLG7\nChLG8\nChLG9\nChLG10\n")
    f.close()
    cmdline = []
    cmdline.append("gatk GenomicsDBImport")
    cmdline.append("--sample-name-map gvcf_map.csv")
    cmdline.append("--genomicsdb-workspace-path consolidated_gvcfs")
    cmdline.append("--intervals interval.list")
    os.system(" ".join(cmdline))

### Export GenomicsDB of consolidated per-sample variants ###
def export_consolidated_database( bucket ):
    os.system("cp -r consolidated_gvcfs/* /scratch/project_2003847/Tribolium_castaneum_SNP_call/consolidated_gvcfs/.")

### Import GenomicsDB of consolidated per-sample variants ###
def import_consolidated_database( bucket ):
    os.system("mkdir consolidated_gvcfs")
    os.system("cp -r /scratch/project_2003847/Tribolium_castaneum_SNP_call/consolidated_gvcfs .")

### Run GATK GenotypeGVCFs ###
def run_gatk_GenotypeGVCFs( population, version ):
    cmdline = []
    cmdline.append("gatk GenotypeGVCFs")
    cmdline.append("-R Tribolium_castaneum_"+version+".fna")
    cmdline.append("-V gendb://consolidated_gvcfs")
    cmdline.append("-O Tribolium_castaneum_"+population+"_"+version+".vcf")
    cmdline.append("--max-genotype-count 1500")
    os.system(" ".join(cmdline))

### Export GenomicsDB database ###
def export_joint_call( bucket, population, version ):
    os.system("cp Tribolium_castaneum_"+population+"_"+version+".vcf /scratch/project_2003847/Tribolium_castaneum_SNP_call/joint_call/Tribolium_castaneum_"+population+"_"+version+".vcf")
    os.system("cp Tribolium_castaneum_"+population+"_"+version+".vcf.idx /scratch/project_2003847/Tribolium_castaneum_SNP_call/joint_call/Tribolium_castaneum_"+population+"_"+version+".vcf.idx")


##################
#      MAIN      #
##################

if __name__ == '__main__':
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(os.environ['DATADIR'])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Import the list of samples   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Import the list of samples")
    import_sample_list(config["bucket"], config["population"], config["version"])
    samples = load_sample_list(config["population"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Run GATK GenomicsDBImport    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if config["task"] in ["GENOMICS_DB", "BOTH"]:
        ### 3.1) Generate the GVCF sample map ###
        print(">> Generate the GVCF sample map")
        generate_gvcf_sample_map(samples)
        ### 3.2) Import all the per-sample GVCF (and their indices) ###
        print(">> Import the GVCF folder")
        import_gvcf_folder(config["bucket"])
        check_gvcf_files(samples)
        ### 3.3) Import the reference genome and compute indices files ###
        print(">> Import reference genome and create indices")
        import_reference_genome(config["bucket"], config["version"])
        create_reference_genome_indices(config["version"])
        ### 3.4) Run GATK GenomicsDBImport ###
        print(">> Run GATK GenomicsDBImport")
        run_gatk_GenomicsDBImport(config["version"])
        ### 3.5) Export the GenomicsDB database ###
        print(">> Export database")
        export_consolidated_database(config["bucket"])
        ### 3.6) Delete the per-sample GVCFs folder ###
        print(">> Delete GVCFs folder")
        delete_gvcf_folder()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Run GATK GenotypeGVCFs       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if config["task"] in ["GENOTYPE_GVCFS", "BOTH"]:
        ### 4.1) Import GenomicsDB database and reference genome if needed ###
        if config["task"] == "GENOTYPE_GVCFS":
            print(">> Import database")
            import_consolidated_database(config["bucket"])
            print(">> Import reference genome and create indices")
            import_reference_genome(config["bucket"], config["version"])
            create_reference_genome_indices(config["version"])
        ### 4.2) Run GATK GenotypeGVCFs ###
        print(">> Run GATK GenotypeGVCFs")
        run_gatk_GenotypeGVCFs(config["population"], config["version"])
        ### 4.3) Export the GenomicsDB database ###
        print(">> Export joint-call GVCF")
        export_joint_call(config["bucket"], config["population"], config["version"])

    print(">> Done.")

