#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_BamPreProcessing.py
# ---------------------
# Run BAM file pre-processing tasks.
# (HPC SCRIPT --> array wrapper)
#
# 1) Bam pre-processing:
#   1.1) Import the unedited BAM file from Allas,
#   1.2) Import the reference genome and generate indices,
#   1.3) Add the read group,
#   1.4) Uncompress the BAM file to SAM,
#   1.5) Edit the SAM file to recalibrate MAPQ values (255 --> 60),
#   1.6) Compress edited SAM file to BAM,
#   1.7) Copy a version of the BAM file to mark duplicates,
#   1.8) Generate BAI index file,
#   1.9) Handle splicing events,
#   1.10) Export edited BAM files to Allas.
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
    parser.add_argument("--population", "-population", help="Samples list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--task", "-task", choices=['COUNTS', 'GATK', 'BOTH'], help="Which task(s) to run")
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

### Import unedited BAM file ###
def import_unedited_bam( input_bam, output_bam, bucket ):
    os.system("a-get "+bucket+"/unedited_bam/"+input_bam)
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

### Sort the BAM file ###
def sort_bam( input_bam, output_bam ):
    os.system("samtools sort "+input_bam+" -o "+output_bam)
    assert_creation(output_bam)
    os.system("rm "+input_bam)
    assert_deletion(input_bam)

### Add the readgroup to the BAM file ###
def add_readgroup( input_bam, output_bam, sample_name ):
    os.system("samtools addreplacerg -r ID:"+sample_name+" -r SM:"+sample_name+" -r PL:ILLUMINA -r LB:Lexogen_mRNA_stranded -o "+output_bam+" "+input_bam)
    assert_creation(output_bam)
    os.system("rm "+input_bam)
    assert_deletion(input_bam)

### Check the presence of the readgroup ###
def check_readgroup( input_bam ):
    proc   = subprocess.Popen(["samtools view -H "+input_bam+" | grep ^\@RG"], stdout=subprocess.PIPE, shell=True)
    output = proc.stdout.read().decode('utf8')
    assert len(output.strip("\n")) > 0, ">> Readgroup has not been written. Exit."
    print(">> Readgroup: "+output.strip("\n"))

### Uncompress BAM file to SAM file ###
def bam_to_sam( input_bam, output_sam ):
    os.system("samtools view -h "+input_bam+" > "+output_sam)
    assert_creation(output_sam)
    os.system("rm "+input_bam)
    assert_deletion(input_bam)

### Replace every MAPQ==255 by 60 ###
def edit_mapq( input_sam, output_sam ):
    f = open(input_sam, "r")
    g = open(output_sam, "w")
    l = f.readline()
    while l:
        if l.startswith("@"):
            g.write(l)
        else:
            split_l = l.strip("\n").split("\t")
            if split_l[4] == "255":
                split_l[4] = "60"
                g.write("\t".join(split_l)+"\n")
            else:
                g.write(l)
        l = f.readline()
    f.close()
    g.close()
    assert_creation(output_sam)
    os.system("rm "+input_sam)
    assert_deletion(input_sam)

### Compress SAM file to BAM file ###
def sam_to_bam( input_sam, output_bam ):
    os.system("samtools view -b "+input_sam+" > "+output_bam)
    assert_creation(output_bam)
    os.system("rm "+input_sam)
    assert_deletion(input_sam)

### Mark duplicates ###
def mark_duplicates( input_bam, output_bam ):
    os.system("picard MarkDuplicates INPUT="+input_bam+" OUTPUT="+output_bam+" METRICS_FILE=metrics.txt USE_JDK_DEFLATER=true USE_JDK_INFLATER=true")
    assert_creation(output_bam)
    os.system("rm "+input_bam)
    assert_deletion(input_bam)

### Handle splicing events ###
def run_SplitNCigarReads( input_bam, output_bam, version ):
    os.system("mkdir tmp_dir")
    reference_genome = "Tribolium_castaneum_"+version+".fna"
    cmdline          = []
    cmdline.append("gatk SplitNCigarReads")
    cmdline.append("-R "+reference_genome)
    cmdline.append("-I "+input_bam)
    cmdline.append("-O "+output_bam)
    cmdline.append("--tmp-dir tmp_dir")
    cmdline.append("-use-jdk-deflater")
    cmdline.append("-use-jdk-inflater")
    os.system(" ".join(cmdline))
    assert_creation(output_bam)
    os.system("rm "+input_bam)
    assert_deletion(input_bam)
    os.system("rm -rf tmp_dir")

### Export BAM file to Allas ###
def export_bam( input_bam, output_bam, bucket ):
    os.system("rclone copyto "+input_bam+" allas:"+bucket+"/bam/"+output_bam)

### Import BAM file ###
def import_bam( input_bam, output_bam, bucket ):
    os.system("a-get "+bucket+"/bam/"+input_bam)
    assert_creation(input_bam)
    if input_bam != output_bam:
        os.system("mv "+input_bam+" "+output_bam)
        assert_deletion(input_bam)
        assert_creation(output_bam)

### Create BAM file index ###
def create_bam_index( input_bam ):
    os.system("samtools index "+input_bam)
    assert_creation(input_bam+".bai")

### Run counts pre-processing tasks ###
def run_counts_tasks( sample_name, config ):
    ### 1) Add the read group ###
    print(">> Add the read group")
    add_readgroup(sample_name+"_unedited.bam", sample_name+"_RG.bam", sample_name)
    check_readgroup(sample_name+"_RG.bam")
    ### 2) Uncompress BAM to SAM ###
    print(">> Uncompress BAM to SAM")
    bam_to_sam(sample_name+"_RG.bam", sample_name+"_RG.sam")
    ### 3) Edit mapping scores ###
    print(">> Edit mapping scores")
    edit_mapq(sample_name+"_RG.sam", sample_name+"_recalibrated.sam")
    ### 4) Compress SAM to BAM ###
    print(">> Compress SAM to BAM")
    sam_to_bam(sample_name+"_recalibrated.sam", sample_name+"_recalibrated.bam")
    ### 5) Create a copy of the BAM file ###
    print(">> Create a copy of the BAM file")
    os.system("cp "+sample_name+"_recalibrated.bam "+sample_name+"_FC.bam")
    ### 6) Export the edited BAM file ###
    print(">> Export edited BAM file")
    export_bam(sample_name+"_FC.bam", sample_name+"_FC.bam", config["bucket"])

### Run GATK pre-processing tasks ###
def run_GATK_tasks( sample_name, config ):
    ### 1) Add the read group ###
    print(">> Add the read group")
    add_readgroup(sample_name+"_unedited.bam", sample_name+"_RG.bam", sample_name)
    check_readgroup(sample_name+"_RG.bam")
    ### 2) Uncompress BAM to SAM ###
    print(">> Uncompress BAM to SAM")
    bam_to_sam(sample_name+"_RG.bam", sample_name+"_RG.sam")
    ### 3) Edit mapping scores ###
    print(">> Edit mapping scores")
    edit_mapq(sample_name+"_RG.sam", sample_name+"_recalibrated.sam")
    ### 4) Compress SAM to BAM ###
    print(">> Compress SAM to BAM")
    sam_to_bam(sample_name+"_recalibrated.sam", sample_name+"_recalibrated.bam")
    ### 5) Create BAM index ###
    print(">> Create BAM index")
    create_bam_index(sample_name+"_recalibrated.bam")
    ### 6) Mark read duplicates on first BAM ###
    print(">> Mark duplicates")
    mark_duplicates(sample_name+"_recalibrated.bam", sample_name+"_duplicates.bam")
    os.system("rm "+sample_name+"_recalibrated.bam.bai")
    ### 7) Create BAM index ###
    print(">> Create BAM index")
    create_bam_index(sample_name+"_duplicates.bam")
    ### 8) Handle splicing events ###
    print(">> Handle splicing events")
    run_SplitNCigarReads(sample_name+"_duplicates.bam", sample_name+"_cigar.bam", config["version"])
    os.system("rm "+sample_name+"_duplicates.bam.bai")
    ### 9) Export the edited BAM file ###
    print(">> Export edited BAM file")
    export_bam(sample_name+"_cigar.bam", sample_name+".bam", config["bucket"])


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
    # 3) Import files from Allas                       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    ### 3.1) Import the unedited BAM file ###
    print(">> Import the unedited BAM file")
    sample_name = sample["sample"]
    import_unedited_bam(sample_name+".bam", sample_name+"_unedited.bam", config["bucket"])

    ### 3.2) Import the reference genome and create index ###
    print(">> Import reference genome and create indices")
    import_reference_genome(config["bucket"], config["version"])
    create_reference_genome_indices(config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Run BAM pre-processing                        #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if config["task"] in ["COUNTS", "BOTH"]:
        run_counts_tasks(sample_name, config)

    elif config["task"] in ["GATK", "BOTH"]:
        run_GATK_tasks(sample_name, config)

    print(">> Done.")

