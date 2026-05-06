from Bio import SeqIO
import datetime as dt
import csv
import os
from collections import defaultdict
import sys
import error_checks as error_checks

def add_bools_to_config(config, fixed_tree, empirical_tree, starting_tree, glm, epoch, verbose):

    if empirical_tree:
        config["fixed_tree"] = True
    else:
        config["fixed_tree"] = fixed_tree
    
    config["empirical_tree"] = empirical_tree
    config["starting_tree"] = starting_tree
    config["glm"] = glm
    config["epoch"] = epoch
    config["verbose"] = verbose

    return config

def decimal_date(date_string):

    if len(date_string.split("-")) == 3:
        date = dt.datetime.strptime(date_string, "%Y-%m-%d").date()
        uncertainty = 0.0
    elif len(date_string.split("-")) == 2:
        date = dt.datetime.strptime(date_string, "%Y-%m").date()
        uncertainty = 0.08333333333333333
    elif len(date_string.split("-")) == 1:
        date = dt.datetime.strptime(date_string, "%Y").date()
        uncertainty = 1.0
    else:
        sys.stderr("date not in recognised format")
        sys.exit(-1)
    
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    final_date = date.year + float(date.toordinal() - start) / year_length

    return str(final_date), uncertainty


def parse_fasta(fasta_list, codon_partitioning):

    fastas = fasta_list.split(",")
    if codon_partitioning:
        cps = codon_partitioning.split(",")
    else:
        cps = ""

    taxa = set()
    if len(fastas) > 1:
        for fasta in fastas:
            for seq in SeqIO.parse(fasta,"fasta"):
                taxa.add(seq.id)
    else:
        for seq in SeqIO.parse(fastas[0],"fasta"):
            taxa.add(seq.id)

    fasta_info = defaultdict(dict)
    if cps != "":
        for fasta,cp in zip(fastas, cps):
            name = fasta.split("/")[-1].strip(".fasta")
            fasta_info[name]["sequences"] = pull_sequences(fasta)
            if cp == "1":
                fasta_info[name]["codon_partitioning"] = True 
            else:
                fasta_info[name]["codon_partitioning"] = False     
    else:
        for fasta in fastas:
            name = fasta.split("/")[-1].strip(".fasta")
            fasta_info[name]["sequences"] = pull_sequences(fasta)
            fasta_info[name]["codon_partitioning"] = False
    
    for name, info in fasta_info.items():
        lst = []
        for i in info["sequences"]:
            lst.append(i.id)

        fasta_info[name]["taxon_list"] = lst
    
    return fasta_info

def pull_sequences(fasta):

    seq_list = []
    for seq in SeqIO.parse(fasta, 'fasta'):
        seq_list.append(seq)

    return seq_list

    
          







        

    




