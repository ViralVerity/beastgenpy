import sys
import csv
import os
import datetime as dt

def check_file_exists(file):

    if not os.path.exists(file):
        sys.stderr.write(f"can't find {file}\n")
        sys.exit(-1)

    return

def check_dates_in_names(config):

    for analysis, info in config["sequence_info"].items():
        for name in info["taxon_list"].items():
            if not "|" in name:
                sys.stderr(f"sequence name format needs to be '|[date]'. {name} in fasta/tree {name} is incorrect \n")
                sys.exit(-1)
            else:
                date_string = name.split("|")[-1]
                if len(date_string.split("-")) == 3:
                    try:
                        date = dt.datetime.strptime(date_string, "%Y-%m-%d").date()
                    except:
                        sys.stderr(f"date in {name} in analysis {analysis} not in recognised format")
                        sys.exit(-1)
                elif len(date_string.split("-")) == 2:
                    try:
                        date = dt.datetime.strptime(date_string, "%Y-%m").date()
                    except:
                        sys.stderr(f"date in {name} in analysis {analysis} not in recognised format")
                        sys.exit(-1)
                elif len(date_string.split("-")) == 1:
                    try:
                        date = dt.datetime.strptime(date_string, "%Y").date()
                    except:
                        sys.stderr(f"date in {name} in analysis {analysis} not in recognised format")
                        sys.exit(-1)
                else:
                    sys.stderr(f"date in {name} in analysis {analysis} not in recognised format")
                    sys.exit(-1)
    return

def check_names_models(config, thisdir):

    path_to_templates = os.path.join(thisdir, "templates")

    allowed_models = set()

    for direc in ("population_models", "substitution_models", "clock_models", "phylogeog_components"):
        for file in os.listdir(os.path.join(path_to_templates, direc)):
            if file.endswith(".xml")
                model = file.split(".")[0]
                allowed_models.add(model)

    for key in ["phylogeography", "subs_model", "clock_model", "population_model"]:
        if config[key] not in allowed_models:
            sys.stderr.write(f"{config[key]} not allowed (either typo or not yet implemented). See docs for implemented models")
            sys.exit(-1)

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

def check_seqs_present(config):

    for fasta_info in config["fasta"]:
        for seq in fasta_info["sequences"]:
            if seq.id not in config["trait_dict"]:
                sys.stderr.write(f"{seq.id} not in trait file\n")
                sys.exit(-1)

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
