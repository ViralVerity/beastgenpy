import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import argparse
import os
import csv
from mako.template import Template
from collections import defaultdict
import glm_funcs as glm_funcs

parser = argparse.ArgumentParser(add_help=False, description='make XML')

parser.add_argument("--fasta")
parser.add_argument("--id-file", dest="id_file")
parser.add_argument("--template")
parser.add_argument("--map-file", dest="map_file")
parser.add_argument("--traits")
# parser.add_argument("--trait-order", dest="trait_order")
#parser.add_argument("--codon-partitioning", dest="codon_partitioning")
parser.add_argument("--trait-file", dest="trait_file")
parser.add_argument("--fixed-tree-file", dest="fixed_tree_file") 
parser.add_argument("--fixed-tree-dir", dest='fixed_tree_dir')
parser.add_argument("--glm", action="store_true")
parser.add_argument("--predictor-info-file", dest="predictor_info_file") #list of predictors that you want logged and standardise
parser.add_argument("--asymmetric-predictor-file", dest="asymmetric_predictor_file")
parser.add_argument("--predictors-dir", dest="predictors_dir")
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
glm = args.glm
predictor_info_file = args.predictor_info_file
directional_predictor_file = args.asymmetric_predictor_file
predictors_dir = args.predictors_dir

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

##move all these bits to core functions or glm functions at some point

#######TAXA IF NO FASTA##############
id_list = []
if id_file:
    with open(id_file) as f:
        for l in f:
            id_list.append(l.strip("\n"))

#######DISCRETE TRAITS#########
new_traits = []
trait_locs = {}
all_trait_options_prep = defaultdict(set)
all_trait_options = defaultdict(list)
if traits:
    new_traits = traits.split(",")
    for number, trait in enumerate(new_traits):
        trait_locs[trait] = int(number)
    
    if trait_file: #ie row per sequence
        trait_dict = defaultdict(list)
        with open(trait_file) as f:
            data = csv.DictReader(f)
            for line in data:
                name = line["sequence_name"]
                for trait in new_traits:
                    trait_dict[name].append(line[trait])
                    all_trait_options_prep[trait].add(line[trait])
                
    for trait, options in all_trait_options_prep.items():
        new_lst = sorted(options)
        all_trait_options[trait] = new_lst
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
    tree_name = ""

################################

#####GLM section######
if glm:
    #need to make this robust to working directory
    re_matrices = glm_funcs.make_twoway_REmatrices(all_trait_options)
    bin_probs = glm_funcs.calculate_binomial_likelihood(all_trait_options)
    trait_to_predictor = defaultdict(dict)
    if len(new_traits) == 1:
        trait_name = new_traits[0]
        actual_predictor_dir = predictors_dir
        trait_to_predictor = glm_funcs.loop_for_processing(actual_predictor_dir, predictor_info_file, directional_predictor_file, trait_name, trait_to_predictor, all_trait_options)

    else:
        for trait_name in new_traits:
            actual_predictor_dir = os.path.join(predictors_dir,trait)
            trait_to_predictor = glm_funcs.loop_for_processing(actual_predictor_dir, predictor_info_file, directional_predictor_file, trait_name, trait_to_predictor, all_trait_options)



mytemplate = Template(filename=template, strict_undefined=True)
f = open(f"{file_stem}.xml", 'w')
f.write(mytemplate.render(id_list=id_list, tree=fixed_tree, tree_name=tree_name, tree_dict=tree_dict, traits=new_traits, trait_dict=trait_dict, trait_locs=trait_locs, all_trait_options=all_trait_options, trait_to_predictor=trait_to_predictor, bin_probs=bin_probs, re_matrices=re_matrices, chain_length=chain_length, log_every=log_every, file_stem=file_stem))
f.close()

