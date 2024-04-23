#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_CreateBamMap.py
# -----------------
# Collect and list BAM files.
# (LOCAL SCRIPT)
#***************************************************************************

import os
import sys
import csv
import argparse
import subprocess

### Parse command line arguments ###
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-path", "-repository-path", help="Repository path")
    parser.add_argument("--source-path", "-source-path", help="BAM files original source path")
    args = parser.parse_args()
    return(vars(args))

### Load the family table ###
def load_family_table( filename ):
    family_table = {}
    file         = open(filename, "r")
    csvreader    = csv.reader(file, delimiter=";")
    header       = next(csvreader)
    for row in csvreader:
        fem               = row[0]
        line              = row[1]
        sire              = row[2]
        dam               = row[3]
        family_table[fem] = {"line":line, "sire":sire, "dam":dam.strip(".A")}
    file.close()
    return family_table

### Load the sequencing runs ###
def load_sequencing_runs( filename ):
    index           = 1
    key_index       = {}
    sequencing_runs = {}
    f               = open(filename, "r")
    l               = f.readline()
    while l:
        date   = l.strip("\n").split(".")[0]
        sample = l.strip("\n").split(".")[1].strip("A-").split("_R1")[0]
        if sample == "L6CT-8-1_HD":
            sample = "L6CT-8-1HD"
        if date not in key_index.keys():
            key_index[date] = str(index)
            index += 1
        assert sample not in sequencing_runs.keys()
        sequencing_runs[sample] = {"date":date, "index":key_index[date]}
        l = f.readline()
    f.close()
    return sequencing_runs

### Load the fitness data from G1 and G21 ###
def load_fitness_data( filename_G1, filename_G21 ):
    fitness_data = {}
    #--------------------------#
    # 1) Load G1 fitness data  #
    #--------------------------#
    file      = open(filename_G1, "r")
    csvreader = csv.reader(file, delimiter="\t")
    header    = next(csvreader)
    for row in csvreader:
        sample_name               = row[0].split("-")[0]+"-"+row[0].split("-")[1][1]+row[0].split("_")[1]
        fitness                   = row[2]
        fitness_data[sample_name] = {"fitness":fitness, "family":""}
    file.close()
    #--------------------------#
    # 2) Load G21 fitness data #
    #--------------------------#
    mapper    = {}
    file      = open(filename_G21, "r")
    csvreader = csv.reader(file, delimiter="\t")
    header    = next(csvreader)
    for row in csvreader:
        left_side   = row[0].split("-")[0]+row[0].split("-")[1]+"-"+row[4].split("_")[1]
        target_env  = row[3]
        mapper_key  = left_side+"-"+target_env
        if mapper_key not in mapper:
            mapper[mapper_key] = 1
        else:
            mapper[mapper_key] += 1
        sample_name               = row[0].split("-")[0]+row[0].split("-")[1]+"-"+row[4].split("_")[1]+"-"+str(mapper[mapper_key])+target_env
        family                    = row[4].split("_")[1]
        fitness                   = row[5]
        fitness_data[sample_name] = {"fitness":fitness, "family":family}
    file.close()
    key_errors = {
        "L5D-3-3D":"L5-D-3-3D",
        "L2H-8-4H":"L2H8-4H",
        "L2D-3-1H": "L2D3-1H",
        "L5CT-2-2D": "L5-CT-2-2D",
        "L5CT-2-4HD": "L5CT2-4HD",
        "L1HD-3-1CT": "L1HD-3-1-CT",
        "Mx1CT-1-2CT": "Mx1-CT-1-2CT",
        "L3CT-5-4H": "L3-CT-5-4H"
    }
    for item in key_errors.items():
        fitness_data[item[1]] = fitness_data[item[0]]
        del fitness_data[item[0]]
    return fitness_data

### Load experimental batches ###
def load_batches( filename ):
    index     = 1
    key_index = {}
    batches   = {}
    f         = open(filename, "r")
    l         = f.readline()
    l         = f.readline()
    while l:
        l      = l.strip("\n").split("\t")
        p1     = l[0].split("-")[0]
        p2     = l[0].split("-")[1].split("_")[0]
        p3     = l[0].split("-")[1].split("_")[1]
        sample = p1+"-"+p2.strip("f")+p3
        batch  = l[5]
        if batch not in key_index.keys():
            key_index[batch] = str(index)
            index += 1
        assert sample not in batches.keys()
        batches[sample] = {"batch":batch, "index":key_index[batch]}
        l = f.readline()
    f.close()
    return batches


##################
#      MAIN      #
##################

if __name__ == '__main__':
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments                        #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(config["repository_path"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Load family table, sequencing runs and fitness data #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load family table, sequencing runs and fitness data")
    FAMILY_TABLE    = load_family_table("data/experiment_data/FamilyTable.csv")
    SEQUENCING_RUNS = load_sequencing_runs("data/experiment_data/fastq_list.txt")
    FITNESS_DATA    = load_fitness_data("data/experiment_data/fitness_data_G1.txt", "data/experiment_data/fitness_data_G21.txt")
    BATCHES         = load_batches("data/experiment_data/Fitness.txt")
    del FAMILY_TABLE['']
    FAMILY_TABLE['-'] = {'line': '-', 'sire': '-', 'dam': '-'}

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Load the list of individual samples                 #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load the list of individual samples")
    hand_annotation = {"L5CT2-4HD":["CT","HD"], "L2D3-1H":["D","H"], "L2H8-4H":["H","H"]}
    samples         = {}
    sample_count    = {}
    proc            = subprocess.Popen(["ls "+config["source_path"]+" | grep '^STAR\|^Map_STAR' | grep -v MappingTest"], stdout=subprocess.PIPE, shell=True)
    folders         = proc.stdout.read().decode('utf8').strip("\n").split("\n")
    for folder in folders:
        path = config["source_path"]+"/"+folder
        if os.path.isfile(path+"/dataset.tsv"):
            file      = open(path+"/dataset.tsv", "r")
            csvreader = csv.reader(file, delimiter="\t")
            header    = next(csvreader)
            for row in csvreader:
                assert "Tcas_3.0" in row[5] or "Tcas3.30" in row[5] or "Tcas5.2" in row[5]
                sample = row[0]
                if sample not in sample_count.keys():
                    sample_count[sample] = 1
                else:
                    sample_count[sample] += 1
                generation     = "NA"
                fem            = "NA"
                line           = "NA"
                fullsib_family = "NA"
                halfsib_family = "NA"
                source_env     = "NA"
                target_env     = "NA"
                fitness        = FITNESS_DATA[sample]["fitness"]
                run_date       = SEQUENCING_RUNS[sample]["date"]
                run_index      = SEQUENCING_RUNS[sample]["index"]
                bam_path       = path+"/"+sample+".bam"
                bai_path       = path+"/"+sample+".bam.bai"
                ########################################################
                # 3.1) If the individual has the G1 nomenclature       #
                ########################################################
                if len(sample.split("-")) == 2 and sample.split("-")[0].isnumeric():
                    generation     = "1"
                    fem            = sample.split("-")[0]
                    line           = FAMILY_TABLE[fem]["line"]
                    fullsib_family = FAMILY_TABLE[fem]["dam"]
                    halfsib_family = FAMILY_TABLE[fem]["sire"]
                    batch          = BATCHES[sample]["batch"]
                    batch_index    = BATCHES[sample]["index"]
                    t1             = sample.split("-")[1]
                    if sample in hand_annotation.keys():
                        source_env = hand_annotation[sample][0]
                        target_env = hand_annotation[sample][1]
                    elif t1.endswith("CT"):
                        source_env = "CT"
                        target_env = "CT"
                    elif t1.endswith("H") and not t1.endswith("HD"):
                        source_env = "H"
                        target_env = "H"
                    elif t1.endswith("D") and not t1.endswith("HD"):
                        source_env = "D"
                        target_env = "D"
                    elif t1.endswith("HD"):
                        source_env = "HD"
                        target_env = "HD"
                ########################################################
                # 3.2) Else if the individual has the G21 nomenclature #
                ########################################################
                else:
                    generation     = "21"
                    fullsib_family = FITNESS_DATA[sample]["family"]
                    batch          = "NA"
                    batch_index    = "NA"
                    ### Manage line ###
                    if sample.startswith("L"):
                        line = sample[0:2]
                    elif sample.startswith("Mx"):
                        line = sample[0:3]
                    ### Manage source environment ###
                    t1 = sample.split("-")[0]
                    t2 = sample.split("-")[1]
                    if sample in hand_annotation.keys():
                        source_env = hand_annotation[sample][0]
                    elif (t1.endswith("CT") or t2.startswith("CT")):
                        source_env = "CT"
                    elif (t1.endswith("H") or t2.startswith("H")) and not t1.endswith("HD") and not t2.startswith("HD"):
                        source_env = "H"
                    elif (t1.endswith("D") or t2.startswith("D")) and not t1.endswith("HD") and not t2.startswith("HD"):
                        source_env = "D"
                    elif (t1.endswith("HD") or t2.startswith("HD")):
                        source_env = "HD"
                    ### Manage target environment ###
                    t1 = sample.split("-")[-1]
                    if sample in hand_annotation.keys():
                        target_env = hand_annotation[sample][1]
                    elif t1.endswith("CT"):
                        target_env = "CT"
                    elif t1.endswith("H") and not t1.endswith("HD"):
                        target_env = "H"
                    elif t1.endswith("D") and not t1.endswith("HD"):
                        target_env = "D"
                    elif t1.endswith("HD"):
                        target_env = "HD"
                ########################################################
                # 3.3) Manage annotation if the BAM file exists        #
                ########################################################
                if os.path.isfile(bam_path):
                    ### Manage annotation ###
                    annotation = row[5].split("/")[4]
                    version    = ""
                    fna_path   = ""
                    if "Tcas_3.0" in row[5]:
                        version  = "Tcas3.0"
                        fna_path = "data/reference_genome/Tribolium_castaneum_Tcas3.0/Tribolium_castaneum_Tcas3.0.fna"
                    elif "Tcas3.30" in row[5]:
                        version  = "Tcas3.30"
                        fna_path = "data/reference_genome/Tribolium_castaneum_Tcas3.30/Tribolium_castaneum_Tcas3.30.fna"
                    elif "Tcas5.2" in row[5]:
                        version  = "Tcas5.2"
                        fna_path = "data/reference_genome/Tribolium_castaneum_Tcas3.30/Tribolium_castaneum_Tcas5.2.fna"
                    ### Save the sample ###
                    samples[sample+"|"+str(sample_count[sample])] = {
                        "sample":sample,
                        "generation":generation,
                        "fem":fem,
                        "line":line,
                        "fullsib_family":fullsib_family,
                        "halfsib_family":halfsib_family,
                        "source_env":source_env,
                        "target_env":target_env,
                        "fitness":fitness,
                        "run_date":run_date,
                        "run_index":run_index,
                        "batch":batch,
                        "batch_index":batch_index,
                        "bam_path":bam_path,
                        "bai_path":bai_path,
                        "annotation":annotation,
                        "version":version,
                        "fna_path":fna_path
                    }
                else:
                    #print("> The BAM file "+bam_path+" does not exist")
                    samples[sample+"|"+str(sample_count[sample])] = {
                        "sample":sample,
                        "generation":"NA",
                        "fem":"NA",
                        "line":"NA",
                        "fullsib_family":"NA",
                        "halfsib_family":"NA",
                        "source_env":"NA",
                        "target_env":"NA",
                        "fitness":"NA",
                        "run_date":"NA",
                        "run_index":"NA",
                        "batch":"NA",
                        "batch_index":"NA",
                        "bam_path":"NA",
                        "bai_path":"NA",
                        "annotation":"NA",
                        "version":"NA",
                        "fna_path":"NA"
                    }
            file.close()
        else:
            print("> No dataset in folder "+path)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Save the BAM map                                    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Save the BAM map")
    header  = ""
    header += "sample;"
    header += "generation;"
    header += "fem;"
    header += "line;"
    header += "fullsib_family;"
    header += "halfsib_family;"
    header += "source_env;"
    header += "target_env;"
    header += "fitness;"
    header += "run_date;"
    header += "run_index;"
    header += "batch;"
    header += "batch_index;"
    header += "bam_path;"
    header += "bai_path;"
    header += "annotation;"
    header += "version;"
    header += "fna_path"
    f = open("data/tribolium_bam/bam_map.csv", "w")
    f.write(header+"\n")
    for item in samples.items():
        line  = ""
        line += item[1]["sample"]+";"
        line += item[1]["generation"]+";"
        line += item[1]["fem"]+";"
        line += item[1]["line"]+";"
        line += item[1]["fullsib_family"]+";"
        line += item[1]["halfsib_family"]+";"
        line += item[1]["source_env"]+";"
        line += item[1]["target_env"]+";"
        line += item[1]["fitness"]+";"
        line += item[1]["run_date"]+";"
        line += item[1]["run_index"]+";"
        line += item[1]["batch"]+";"
        line += item[1]["batch_index"]+";"
        line += item[1]["bam_path"]+";"
        line += item[1]["bai_path"]+";"
        line += item[1]["annotation"]+";"
        line += item[1]["version"]+";"
        line += item[1]["fna_path"]
        f.write(line+"\n")
    f.close()

    print(">> Done.")

