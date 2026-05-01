import sys
import csv
import os

def check_seqs_present(config):

    for fasta_info in config["fasta"]:
        for seq in fasta_info["sequences"]:
            if seq.id not in config["trait_dict"]:
                sys.stderr.write(f"{seq.id} not in trait file\n")
                sys.exit(-1)

def check_alignment(config):
    
#if not fixed tree, check alignment is present
    return


def check_gp_cutoff(config):

    #if skygrid, check gridpoints and cutoff are present

    return


def check_tree_dict(config):

    #if fixed/empirical check tree files are there

    return


def check_names_models(config):

    #that I've implemented skygrid/relaxed clock/whatever
    #matches to existing template blocks

    return

def check_polygon_dir(config):

    #polygon dir is findable 

    return

