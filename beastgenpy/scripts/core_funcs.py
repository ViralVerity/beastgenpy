from Bio import SeqIO
import datetime as dt
import csv
import os
from collections import defaultdict

def add_bools_to_config(config, multi_tree, fixed_tree, dta, glm, continuous_phylogeog):

    config["multi_tree"] = multi_tree
    config["fixed_tree"] = fixed_tree
    config["dta"] = dta
    config["glm"] = glm
    config["continuous_phylogeog"] = continuous_phylogeog

    return config

def decimal_date(date_string):

    date = dt.datetime.strptime(date_string, "%Y-%m-%d").date()
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    final_date = date.year + float(date.toordinal() - start) / year_length

    return str(final_date)


def parse_fasta(fasta_list, codon_partitioning):

    fastas = fasta_list.split(",")
    cps = codon_partitioning.split(",")

    taxa = []
    for seq in SeqIO.parse(fastas[0],"fasta"):
        taxa.append(seq.id)

    fasta_info = []
    for fasta,cp in zip(fastas, cps):
        inner_dict = {}
        inner_dict["name"] = fasta.split("/")[-1].strip(".fasta")

        inner_dict["sequences"] = pull_sequences(fasta)
        if cp == "1":
            inner_dict["codon_partitioning"] = True 
        else:
            inner_dict["codon_partitioning"] = False
        
        fasta_info.append(inner_dict)
 
    return taxa, fasta_info

def pull_sequences(fasta):

    seq_list = []
    for seq in SeqIO.parse(fasta, 'fasta'):
        seq_list.append(seq)

    return seq_list

    
def get_taxa_no_fasta(id_file, id_file_dir, fixed_tree_file, config):
#taxa data structure needs thinking about - mostly the key

    taxa = defaultdict(list)
    if id_file:
        if config["fixed_tree"]:
            name = fixed_tree_file.split("/")[-1].split(".")[0]
        else:
            name = "taxa"
        with open(id_file) as f:
            data = csv.DictReader(f,delimiter="\t")
            for l in data:
                taxa[name].append(l["name"])

    elif id_file_dir: #for taxon sets or for multitree
        for f in os.listdir(id_file_dir):
            if f.endswith(".csv"):
                file_name = f.split("/")[-1].split(".")[0]
                with open(os.path.join(id_file_dir,f)) as open_f:
                    data = csv.DictReader(open_f)
                    for l in data:
                        taxa[file_name].append(l["name"])
    else:
        taxa = False #not sure if I want this - shouldn't it just error out? there's now way to get IDs if not those args or fasta

    return taxa

#some kind of taxon set function

            
def fixed_tree_parsing(fixed_tree_file, fixed_tree_dir, config):
#errors:
#check dir is there if multitree
#need one of these args if fixed tree is selected

    tree_file_dict = {}
    newick_dict = {}
    if config["multi_tree"]:
        for f in os.listdir(fixed_tree_dir):
            if f.endswith(".newick") or f.endswith(".nwk") or f.endswith("tree"): #for now, will only take newick strings
                tree_name = f.split("/")[-1].split(".")[0]
                with open(os.path.join(fixed_tree_dir,f)) as open_f:
                    for l in open_f:
                        fixed_tree = l.strip("\n").lstrip("[&R] ")
                    newick_dict[tree_name] = fixed_tree
                tree_file_dict[tree_name] = f #do I need this data structure?
        tree_name = "" #not sure why this is here

    else:
        with open(fixed_tree_file) as f:
            for l in f:
                fixed_tree = l.strip("\n").lstrip("[&R] ")

        tree_name = fixed_tree_file.split("/")[-1].split(".")[0]
        tree_file_dict[tree_name] = fixed_tree_file.split("/")[-1]
        newick_dict[tree_name] = fixed_tree

    return tree_name, tree_file_dict, newick_dict
    
       







        

    




