import os
from collections import defaultdict
import sys
import dendropy
import error_checks as error_checks

def parse_fixed_trees(config, tree_file, tree_dir):
    
    tree_info = defaultdict(dict)
    if tree_file:
        error_checks.check_tree_file(tree_file)
        tree_info = do_tree_parse(config, tree_file, tree_info)

    elif tree_dir:
        for file in os.listdir(tree_dir):
            tree_file = os.path.join(tree_dir, file)
            error_checks.check_tree_file(tree_file)
            tree_info = do_tree_parse(config, tree_file, tree_info)

    if not config["sequence_info"]:
        config["sequence_info"] = defaultdict(dict)

    for name, info in tree_info.items():
        for key, value in info.items():
            config["sequence_info"][name][key] = value
    
    return config
   
def do_tree_parse(config, tree_file, tree_info):
    
    name = tree_file.split("/")[-1].strip(".tree").strip(".treefile").strip(".trees").strip(".newick")
    tree_info[name]["tree_file"] = tree_file
    if not config["sequence_info"]:
        taxa = get_taxa_from_treefile(config, tree_file)
        tree_info[name]["taxon_list"] = taxa
    
    return tree_info

def get_taxa_from_treefile(config,tree_file):
    
    if config["verbose"]:
        sys.stdout.write("reading trees in to get taxa list, this might take a minute\n")
    
    if tree_file.endswith(".trees") or tree_file.endswith(".nexus"):
        schema = "nexus"
    else:
        schema = "newick"

    tree = dendropy.Tree.get(path=tree_file, schema=schema, default_rooted=True, preserve_underscores=True, offset=0)            
    taxa = []
    for taxon in tree.preorder_leaf_iter():
        taxa.append(taxon.taxon.label)
    return taxa


def parse_starting_trees(config, tree_file, tree_dir):

    tree_dict = {}
    if tree_file:
        error_checks.check_tree_file(tree_file, format="newick")
        name = tree_file.split("/")[-1].strip(".tree").strip(".treefile").strip(".trees").strip(".newick")
        tree = dendropy.Tree.get(path=tree_file, schema="newick", default_rooted=True, preserve_underscores=True)            
        config["sequence_info"][name]["tree_string"] = tree.extract_tree()

    elif tree_dir:
        for file in os.listdir(tree_dir):
            path = os.path.join(tree_dir, file)
            error_checks.check_tree_file(path, format="newick")
            tree = dendropy.Tree.get(path=tree_file, schema="newick", default_rooted=True, preserve_underscores=True)            
            config["sequence_info"][name]["tree_string"] = tree.extract_tree()

    return config