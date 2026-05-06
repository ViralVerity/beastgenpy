import sys
import csv
import os

def check_seqs_present(config):

    for fasta_info in config["fasta"]:
        for seq in fasta_info["sequences"]:
            if seq.id not in config["trait_dict"]:
                sys.stderr.write(f"{seq.id} not in trait file\n")
                sys.exit(-1)

def check_file_exists(file):

    #alignment, trees, trait file

    return

def check_gp_cutoff(config, cutoff, gridpoints):

    if not cutoff:
        sys.stderr.write(f"no cutoff specified for skygrid model\n")
        sys.exit(-1)
    else:
        config["cutoff"] = int(cutoff)
    
    if not gridpoints:
        if config["verbose"]:
            sys.stderr.write(f"no number of gridpoints specified for skygrid model, using one per year\n")
        config["gridpoints"] = int(cutoff)
    
    return config

def check_dates_in_names(config):


    return

def check_names_models(config):

    #that I've implemented skygrid/relaxed clock/whatever
    #matches to existing template blocks

    return

def check_polygon_dir(config):

    #polygon dir is findable 

    return

def check_phylogeog_value(config, phylogeography):

    allowed = ["continuous", "discrete"]
    if phylogeography not in allowed:
        sys.stderr.write(f"{phylogeography} not allowed option, allowed options are discrete or continuous")
        sys.exit(-1)
    else:
        config["phylogeography"] = phylogeography

    return config

def check_headers_trait(config):

    with open(config["trait_file"]):
        data = csv.DictReader(f)
        headers = data.fieldnames

    if "latitude" not in headers:
        sys.stderr.write("latitude missing from headers in continuous trait file\n")
        sys.exit(-1)
        
    if "longitude" not in headers:
        sys.stderr.write("longitude missing from headers in continuous trait file\n")
        sys.exit(-1)

    if "coordinates" not in headers: 
        sys.stderr.write("coordinates missing from headers in continuous trait file\n")
        sys.exit(-1)
