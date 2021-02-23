import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import argparse
import os
from mako.template import Template

parser = argparse.ArgumentParser(add_help=False, description='make XML')

parser.add_argument("--fasta")
parser.add_argument("--map-file", dest="map_file")
parser.add_argument("--traits")
parser.add_argument("--trait-order", dest="trait_order")
parser.add_argument("--codon-partitioning", dest="codon_partitioning")
parser.add_argument("--fixed-tree-file", dest="fixed_tree_file")

args = parser.parse_args()

fasta_file = args.fasta
maps = args.map_file
traits = args.traits
trait_order = args.trait_order
codon_partitioning = args.codon_partitioning
fixed_tree_file = args.fixed_tree_file

new_traits = []
trait_locs = {}
if traits:
    new_traits = traits.split(",")
    for number, trait in zip(trait_order, new_traits):
        trait_locs[trait] = int(number)

if fixed_tree_file:
    with open(fixed_tree_file) as f:
        for l in f:
            fixed_tree = l.strip("\n")

print(fixed_tree)


mytemplate = Template(filename='./phylogeography_template.txt', strict_undefined=True)
f = open("test.xml", 'w')
f.write(mytemplate.render(fasta=fasta_file, traits=new_traits, trait_locs=trait_locs, codon_partitioning=codon_partitioning))
f.close()

