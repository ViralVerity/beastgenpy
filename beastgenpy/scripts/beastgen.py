import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import argparse
import os
import csv
from mako.template import Template
from collections import defaultdict
import glm_funcs as glm_funcs
import core_funcs as core_funcs

parser = argparse.ArgumentParser(add_help=False, description='make XML')

parser.add_argument("--fasta")
parser.add_argument("--id-file", dest="id_file")
parser.add_argument("--id-file-dir", dest="id_file_dir")
parser.add_argument("--template", required=True)
parser.add_argument("--map-file", dest="map_file")
parser.add_argument("--traits")
# parser.add_argument("--trait-order", dest="trait_order")
#parser.add_argument("--codon-partitioning", dest="codon_partitioning")
parser.add_argument("--cont-phylo", action="store_true",dest="cont_phylo")
parser.add_argument("--trait-file", dest="trait_file")
parser.add_argument("--fixed-tree-file", dest="fixed_tree_file") 
parser.add_argument("--fixed-tree-dir", dest='fixed_tree_dir')
parser.add_argument("--glm", action="store_true")
parser.add_argument("--predictor-info-file", dest="predictor_info_file") #list of predictors that you want logged and standardise
parser.add_argument("--asymmetric-predictor-file", dest="asymmetric_predictor_file")
parser.add_argument("--taxon-set-file", dest="taxon_set_file") #must be provided if multi-tree and not fixed tree
parser.add_argument("--predictors-dir", dest="predictors_dir")
parser.add_argument("--chainlen", default="100000000")
parser.add_argument("--log", default="10000")
parser.add_argument("--file-stem", dest="file_stem")

args = parser.parse_args()

template = args.template
fasta_file = args.fasta
id_file = args.id_file
id_file_dir = args.id_file_dir
maps = args.map_file
traits = args.traits
# trait_order = args.trait_order
trait_file = args.trait_file
# codon_partitioning = args.codon_partitioning
cont_phylo=args.cont_phylo
fixed_tree_file = args.fixed_tree_file
fixed_tree_dir = args.fixed_tree_dir
chain_length = args.chainlen
log_every = args.log
file_stem = args.file_stem
glm = args.glm
predictor_info_file = args.predictor_info_file
directional_predictor_file = args.asymmetric_predictor_file
taxon_sets_file = args.taxon_set_file
predictors_dir = args.predictors_dir

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

##move all these bits to core functions or glm functions at some point. Would be nice to have different modules for different analysis types
#realistically, a lot of this needs re-writing and re-framing. Looks like the actual functions are ok, but some of the data structures need thinking about
config = {} #would be nice to set this up like civet eventually - even if not defaults, at least do initialise all the keys
trait_dict = {}

#######TAXA IF NO FASTA##############

taxa = defaultdict(list)
if id_file:
    with open(id_file) as f:
        data = csv.DictReader(f,delimiter="\t")
        for l in data:
            taxa["taxa"].append(l["name"])

elif id_file_dir:
    for f in os.listdir(id_file_dir):
        if f.endswith(".csv"):
            file_name = f.split("/")[-1].split(".")[0]
            with open(os.path.join(id_file_dir,f)) as open_f:
                data = csv.DictReader(open_f)
                for l in data:
                    taxa[file_name].append(l["name"])
else:
    taxa = ""


#######DISCRETE TRAITS#########

trait_locs = {}
options_per_tree = defaultdict(dict) #trait to tree to options
all_trait_options_prep = defaultdict(set)
all_trait_options = defaultdict(list)
if traits:
    traits = traits.split(",")
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



    # if taxon_set_file: #if multitree
    #     options_per_tree = core_funcs.get_options_per_tree(trait_locs, trait_dict, taxon_set_file)
#################################

###########FIXED TREE###############
tree_file_dict = {}
if fixed_tree_file:
    with open(fixed_tree_file) as f:
        for l in f:
            fixed_tree = l.strip("\n").lstrip("[&R] ")

    tree_name = fixed_tree_file.split("/")[-1].split(".")[0]
    tree_file_dict[tree_name] = fixed_tree_file.split("/")[-1]
else:
    tree_name = ""
    fixed_tree = ""

newick_dict = {}

if fixed_tree_dir:
    for f in os.listdir(fixed_tree_dir):
        if f.endswith(".newick") or f.endswith(".nwk") or f.endswith("tree"): #for now, will only take newick strings
            tree_name = f.split("/")[-1].split(".")[0]
            with open(os.path.join(fixed_tree_dir,f)) as open_f:
                for l in open_f:
                    fixed_tree = l.strip("\n").lstrip("[&R] ")
                newick_dict[tree_name] = fixed_tree
            tree_file_dict[tree_name] = f
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
else:
    re_matrices = False
    bin_probs = False
    trait_to_predictor = False



#####CONTINUOUS PHYLOGEOGRAPHY####
if cont_phylo:
    #this template should be generalised to be fixed tree or not fixed, and multi/not multi
    traits = ["longitude", "latitude", "coordinates"]
    trait_dict = core_funcs.process_coordinates(trait_file)
    overall_trait = "coordinates"
else:
    overall_trait = ""
    


###Final set up####
config["taxa"] = taxa

#fixed tree analyses
config["tree"] = fixed_tree
config["tree_name"] = tree_name
config["newick_dict"] = newick_dict
config["tree_file_dict"] = tree_file_dict

##Trait analysis
config["trait_dict"] = trait_dict
##continuous trait analysis
config["overall_trait"] = overall_trait
##discrete trait analysis
config["traits"] = traits
config["trait_locs"] = trait_locs
config["all_trait_options"] = all_trait_options
##DTA GLM
config["trait_to_predictor"] = trait_to_predictor
config["bin_probs"]=bin_probs
config["re_matrices"] =re_matrices

##general options
config["chain_length"]=chain_length
config["log_every"]=log_every
config["file_stem"]=file_stem


mytemplate = Template(filename=template, strict_undefined=True)
f = open(f"{file_stem}.xml", 'w')
f.write(mytemplate.render(config=config))
f.close()

