import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import argparse
import os
import csv
from mako.template import Template
from collections import defaultdict

parser = argparse.ArgumentParser(add_help=False, description='make XML')

parser.add_argument("--fasta")
parser.add_argument("--id-file", dest="id_file")
parser.add_argument("--template")
parser.add_argument("--map-file", dest="map_file")
parser.add_argument("--traits")
# parser.add_argument("--trait-order", dest="trait_order")
#parser.add_argument("--codon-partitioning", dest="codon_partitioning")
parser.add_argument("--trait-file", dest="trait_file")
parser.add_argument("--fixed-tree-file", dest="fixed_tree_file") #will need to make this a dir to deal with multi trees
parser.add_argument("--fixed-tree-dir", dest='fixed_tree_dir')
parser.add_argument("--chainlen", default="100000000")
parser.add_argument("--log", default="10000")
parser.add_argument("--file-stem", dest="file_stem")

args = parser.parse_args()

template = args.template
fasta_file = args.fasta
id_file = args.id_file
maps = args.map_file
traits = args.traits
# trait_order = args.trait_order
trait_file = args.trait_file
# codon_partitioning = args.codon_partitioning
fixed_tree_file = args.fixed_tree_file
fixed_tree_dir = args.fixed_tree_dir
chain_length = args.chainlen
log_every = args.log
file_stem = args.file_stem


##move all these bits to core functions at some point

#######TAXA IF NO FASTA##############
id_list = []
if id_file:
    with open(id_file) as f:
        for l in f:
            id_list.append(l.strip("\n"))

#######DISCRETE TRAITS#########
new_traits = []
trait_locs = {}
all_trait_options = defaultdict(set)
if traits:
    new_traits = traits.split(",")
    for number, trait in enumerate(new_traits):
        trait_locs[trait] = int(number)
    
    if trait_file: 
        trait_dict = defaultdict(list)
        with open(trait_file) as f:
            data = csv.DictReader(f)
            for line in data:
                name = line["sequence_name"]
                for trait in new_traits:
                    trait_dict[name].append(line[trait])
                    all_trait_options[trait].add(line[trait])
#################################

###########FIXED TREE###############
if fixed_tree_file:
    with open(fixed_tree_file) as f:
        for l in f:
            fixed_tree = l.strip("\n").lstrip("[&R] ")

    tree_name = fixed_tree_file.split("/")[-1].split(".")[0]
else:
    tree_name = ""
    fixed_tree = ""

tree_dict = {}
if fixed_tree_dir:
    for f in os.listdir(fixed_tree_dir):
        if f.endswith(".newick") or f.endswith(".nwk") or f.endswith("tree"): #for now, will only take newick strings
            with open(os.path.join(fixed_tree_dir,f)) as open_f:
                for l in open_f:
                    fixed_tree = l.strip("\n").lstrip("[&R] ")
                tree_name = f.split("/")[-1].split(".")[0]
                tree_dict[tree_name] = fixed_tree


################################


mytemplate = Template(filename=template, strict_undefined=True)
f = open(f"{file_stem}.xml", 'w')
f.write(mytemplate.render(id_list=id_list, tree=fixed_tree, tree_name=tree_name, tree_dict=tree_dict, traits=new_traits, trait_dict=trait_dict, trait_locs=trait_locs, all_trait_options=all_trait_options, chain_length=chain_length, log_every=log_every, file_stem=file_stem))
f.close()

