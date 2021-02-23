import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import argparse
import os
import csv
from mako.template import Template
from collections import defaultdict

parser = argparse.ArgumentParser(add_help=False, description='make XML')

parser.add_argument("--fasta")
parser.add_argument("--map-file", dest="map_file")
parser.add_argument("--traits")
# parser.add_argument("--trait-order", dest="trait_order")
#parser.add_argument("--codon-partitioning", dest="codon_partitioning")
parser.add_argument("--trait-file", dest="trait_file")
parser.add_argument("--fixed-tree-file", dest="fixed_tree_file") #will need to make this a dir to deal with multi trees

args = parser.parse_args()

fasta_file = args.fasta
maps = args.map_file
traits = args.traits
# trait_order = args.trait_order
trait_file = args.trait_file
# codon_partitioning = args.codon_partitioning
fixed_tree_file = args.fixed_tree_file


#######DISCRETE TRAITS#########
new_traits = []
trait_locs = {}
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
#################################

###########FIXED TREE###############
if fixed_tree_file:
    with open(fixed_tree_file) as f:
        for l in f:
            fixed_tree = l.strip("\n").lstrip("[&R] ")

tree_name = fixed_tree_file.split("/")[-1].split(".")[0]
################################


mytemplate = Template(filename='/Users/s1743989/Documents/GitHub/beastgenpy/beastgenpy/templates/phylogeography_template.txt', strict_undefined=True)
f = open("test.xml", 'w')
f.write(mytemplate.render(fasta=fasta_file, tree=fixed_tree, tree_name=tree_name, traits=new_traits, trait_dict=trait_dict, trait_locs=trait_locs))
f.close()

