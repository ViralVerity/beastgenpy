import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import argparse
import os
import csv
from mako.template import Template
from collections import defaultdict

import glm_funcs as glm_funcs
import core_funcs as core_funcs
import trait_analysis_funcs as trait_funcs

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False, description='make XML')

    tax_group = parser.add_argument_group('Taxa options')
    tax_group.add_argument("--fasta", help="fasta file containing alignment to analyse")
    tax_group.add_argument("--id-file", dest="id_file", help="File containing taxon IDs in csv format if no fasta file provided")
    tax_group.add_argument("--id-file-dir", dest="id_file_dir", help="Directory containing files with sets of taxon IDs in csv format")

    tree_group = parser.add_argument_group("Tree options")
    tree_group.add_argument("--multi-tree", action="store_true", dest="multi_tree", help="Perform joint analysis on multiple monophyletic trees")
    tree_group.add_argument("--fixed-tree", action="store_true", dest="fixed_tree", help="Perform analysis on fixed tree or trees")
    tree_group.add_argument("--fixed-tree-file", dest="fixed_tree_file", help="File containing single fixed tree in newick format") 
    tree_group.add_argument("--fixed-tree-dir", dest='fixed_tree_dir', help="Directory containing multiple fixed trees in newick format")


    trait_group = parser.add_argument_group("trait_analysis_group")
    trait_group.add_argument("--dta", action="store_true", help="Flag to run a discrete trait analysis")
    trait_group.add_argument("--traits", help="Comma separated list of traits for discrete trait analysis")
    trait_group.add_argument("--discrete-trait-file", dest="discrete_trait_file", help="File containing values for discrete traits, one line per sequence, in csv format")
    trait_group.add_argument("--trait-location-in-name", dest="trait_location_in_name", help="Information to pull trait value from taxon label. Should be provided as trait=position_in_name with different traits separated by commas. e.g.location=5, counting from 1")
    trait_group.add_argument("--trait-delimiter", dest="trait_delimiter", default="|", help="Separater in taxon label to pull trait value from taxon label. Default='|'")

    trait_group.add_argument("--glm", action="store_true", help="Flag to run a generalised linear model on discrete trait predictors")
    trait_group.add_argument("--predictor-info-file", dest="predictor_info_file", help="File containing which predictors should be logged and standardised in csv format. Headers must be 'predictor' and 'logged_transformed_and_standardised'") 
    trait_group.add_argument("--directional-predictor-file", dest="directional_predictor_file", help="File containing directional i.e. one-way predictors for a trait in csv format. Headers must be the name of the trait, and then the predictor names") 
    trait_group.add_argument("--predictors-dir", dest="predictors_dir", help="Directory containing matrices for symmetrical predictors. Each file should have the first row and column as trait values, with the predictor value in the appropriate row/column combination.")

    trait_group.add_argument("--continuous-phylogeog", action="store_true",dest="continuous_phylogeog", help="Flag to run a continuous phylogeographic analysis")
    trait_group.add_argument("--continuous-trait-file", dest="continuous_trait_file", help="File containing coordinate values under headers 'taxon,longitude,latitude' for each sequence for continuous phylogeographic analysis")


    gen_group = parser.add_argument_group('General options')
    gen_group.add_argument("--template", required=True, help="Template for making the XML")
    gen_group.add_argument("--chain-length", default="100000000", help="Number of states to run the MCMC chain for. Default=100m")
    gen_group.add_argument("--log", default="10000", help="How often to log tree and log files. Default=10,000")
    gen_group.add_argument("--file-stem", dest="file_stem", help="File stem for analysis")
    gen_group.add_argument("-h","--help",action="store_true",dest="help")


    args = parser.parse_args()

    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)

    config = {} 

    config = core_funcs.add_bools_to_config(config, args.multi_tree, args.fixed_tree, args.dta, args.glm, args.continuous_phylogeography)

    if args.fasta:
        config["taxa"] = core_funcs.parse_fasta(args.fasta) #needs writing
        config["fasta"] = args.fasta
    else:
        config["taxa"] = core_funcs.get_taxa_no_fasta(args.id_file, args.id_file_dir)
        config["fasta"] = False

    if config["fixed_tree"]:
        config["tree_name"], config["tree_file_dict"], config["newick_dict"] = core_funcs.fixed_tree_parsing(args.fixed_tree_file, args.fixed_tree_dir, config)
    else:
        config["tree_name"] = False
        config["tree_file_dict"] = False
        config["newick_dict"] = False


    if config["dta"]:
        config["traits"],config["trait_index"], config["all_trait_options"], config["trait_dict"], config["options_per_tree"] = trait_funcs.parse_discrete_funcs(args.traits, args.discrete_trait_file,  args.trait_loc_in_name_input, args.trait_delimiter, config)
    else:
        config["trait_index"] = False
        config["all_trait_options"] = False
        config["options_per_tree"] = False
        


    if config["glm"]: #needs dta bool to be true
        config["trait_to_predictor"], config["re_matrices"], config["bin_probs"] = glm_funcs.run_glm_functions(args.predictor_dir_input, args.predictor_info_file, args.directional_predictor_file, config)
    else:
        config["trait_to_predictor"] = False
        config["re_matrices"] = False
        config["bin_probs"] = False


    if config["continuous_phylogeog"]:
        #this template should be generalised to be fixed tree or not fixed, and multi/not multi
        config["traits"], config["trait_dict"], config["overall_trait"] = trait_functions.continuous_phylogeography_processing(args.continuous_trait_file)
    else:
        config["overall_trait"] = False


    if not config["continuous_phylogeog"] or config["dta"]:  
        config["traits"] = False
        config["trait_dict"] = False


    ##general options
    config["chain_length"] = chain_length
    config["log_every"] = log_every
    config["file_stem"] = file_stem


    mytemplate = Template(filename=template, strict_undefined=True)
    f = open(f"{file_stem}.xml", 'w')
    f.write(mytemplate.render(config=config))
    f.close()


if __name__ == '__main__':
    main()



