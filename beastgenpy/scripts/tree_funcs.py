import os
from collections import defaultdict
import sys
import dendropy
import error_checks as error_checks

accetable_endings = ["newick", "nwk", "nexus", "nxs", "tree", "trees", "treefile"]

def parse_fixed_trees(config, tree_file, tree_dir):
    
    tree_info = defaultdict(dict)
    if tree_file:
        error_checks.check_tree_file(tree_file)
        tree_info = do_tree_parse(config, tree_file, tree_info)

    elif tree_dir:
        for file in os.listdir(tree_dir):
            if file.split(".")[-1] in accetable_endings:
                tree_file = os.path.join(tree_dir, file)
                error_checks.check_tree_file(tree_file)
                tree_info = do_tree_parse(config, tree_file, tree_info)

    if "sequence_info" not in config:
        config["sequence_info"] = defaultdict(dict)

    for name, info in tree_info.items():
        for key, value in info.items():
            config["sequence_info"][name][key] = value
    
    return config
   
def do_tree_parse(config, tree_file, tree_info):
    
    name = tree_file.split("/")[-1].split(".")[0]
    tree_info[name]["tree_file"] = tree_file
    if "sequence_info" not in config:
        taxa = get_taxa_from_treefile(config, tree_file)
        tree_info[name]["taxon_list"] = taxa
    
    return tree_info

def get_taxa_from_treefile(config,tree_file):
    
    if config["verbose"]:
        sys.stdout.write("reading trees in to get taxa list, this might take a minute\n")
    
    if tree_file.endswith(".trees") or tree_file.endswith(".nexus"): #add a thing saying that tree or treefile or trees is ambiguous so need --format flag
        scheme = "nexus"
    else:
        scheme = "newick"

    tree = dendropy.Tree.get(path=tree_file, schema=scheme, rooting="force-rooted", preserve_underscores=True, tree_offset=0)            
    taxa = []
    for taxon in tree.leaf_node_iter():
        taxa.append(taxon.taxon.label)
    return taxa


def parse_starting_trees(config, tree_file, tree_dir):

    tree_dict = {}
    if tree_file:
        error_checks.check_tree_file(tree_file, format="newick")
        name = tree_file.split("/")[-1].split(".")[0]
        tree = dendropy.Tree.get(path=tree_file, schema="newick",rooting="force-rooted", preserve_underscores=True)            
        config["sequence_info"][name]["tree_string"] = tree.extract_tree()

    elif tree_dir:
        for file in os.listdir(tree_dir):
            if file.split(".")[-1] in accetable_endings:
                path = os.path.join(tree_dir, file)
                error_checks.check_tree_file(path, format="newick")
                tree = dendropy.Tree.get(path=tree_file, schema="newick",rooting="force-rooted", preserve_underscores=True)            
                config["sequence_info"][name]["tree_string"] = tree.extract_tree()

    return config