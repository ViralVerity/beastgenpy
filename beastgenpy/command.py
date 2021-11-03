import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import argparse
import os
import csv
from mako.template import Template
from mako.runtime import Context
from mako.exceptions import RichTraceback
from io import StringIO
from collections import defaultdict

import glm_funcs as glm_funcs
import core_funcs as core_funcs
import trait_analysis_funcs as trait_funcs

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False, description='make XML')

    tax_group = parser.add_argument_group('Taxa options')
    tax_group.add_argument("--fasta", help="comma separated list of fasta files containing alignment to analyse")
    tax_group.add_argument("--codon-partitioning", "-cp", dest="codon_partitioning", help="comma separated list of 1s and 0s for which alignments have codon partitioning")
    tax_group.add_argument("--id-file", dest="id_file", help="File containing taxon IDs in csv format if no fasta file provided")
    tax_group.add_argument("--id-file-dir", dest="id_file_dir", help="Directory containing files with sets of taxon IDs in csv format")

    tree_group = parser.add_argument_group("Tree options")
    tree_group.add_argument("--multi-tree", action="store_true", dest="multi_tree", help="Perform joint analysis on multiple monophyletic trees")
    tree_group.add_argument("--fixed-tree", action="store_true", dest="fixed_tree", help="Perform analysis on fixed tree or trees")
    tree_group.add_argument("--fixed-tree-file", dest="fixed_tree_file", help="File containing single fixed tree in newick format") 
    tree_group.add_argument("--fixed-tree-dir", dest='fixed_tree_dir', help="Directory containing multiple fixed trees in newick format")
    tree_group.add_argument("--starting-tree", action="store_true", dest="starting_tree", help="flag for adding a starting tree")
    tree_group.add_argument("--starting-tree-file", dest='starting_tree_file', help="file containing newick string for starting tree")

    growth_model_group = parser.add_argument_group("Growth model analysis")
    growth_model_group.add_argument("--growth-model", dest="growth_model", default="skygrid")
    growth_model_group.add_argument("--sg-cutoff", dest="sg_cutoff")
    growth_model_group.add_argument("--sg-gridpoints", dest="sg_gridpoints")

    trait_group = parser.add_argument_group("trait_analysis_group")
    trait_group.add_argument("--dta", action="store_true", help="Flag to run a discrete trait analysis")
    trait_group.add_argument("--traits", help="Comma separated list of traits for discrete trait analysis")
    trait_group.add_argument("--discrete-trait-file", dest="discrete_trait_file", help="File containing values for discrete traits, one line per sequence, in csv format")
    trait_group.add_argument("--trait-location-in-name", dest="trait_location_in_name", help="Information to pull trait value from taxon label. Should be provided as trait=position_in_name with different traits separated by commas. e.g.location=5, counting from 1")
    trait_group.add_argument("--trait-delimiter", dest="trait_delimiter", default="|", help="Separater in taxon label to pull trait value from taxon label. Default='|'")
    trait_group.add_argument("--ambiguities", help="tsv containing ambiguities for dta with two columns: one with header 'ambiguity' and one with 'options' containing a comma separated list")

    trait_group.add_argument("--glm", action="store_true", help="Flag to run a generalised linear model on discrete trait predictors")
    trait_group.add_argument("--predictor-info-file", dest="predictor_info_file", help="File containing which predictors should be logged and standardised in csv format. Headers must be 'predictor' and 'logged_transformed_and_standardised'") 
    trait_group.add_argument("--asymmetric-predictor-file", dest="asymmetric_predictor_file", help="File containing asymmetric i.e. one-way predictors for a trait in csv format. Headers must be the name of the trait, and then the predictor names") 
    trait_group.add_argument("--symmetric-predictor-dir", dest="symmetric_predictor_dir", help="Directory containing matrices for symmetrical predictors. Each file should have the first row and column as trait values, with the predictor value in the appropriate row/column combination.")

    trait_group.add_argument("--continuous-phylogeog", action="store_true",dest="continuous_phylogeog", help="Flag to run a continuous phylogeographic analysis")
    trait_group.add_argument("--continuous-trait-file", dest="continuous_trait_file", help="File containing coordinate values under headers 'taxon,longitude,latitude' for each sequence for continuous phylogeographic analysis")


    gen_group = parser.add_argument_group('General options')
    gen_group.add_argument("--template", required=True, help="Template for making the XML")
    gen_group.add_argument("--chain-length", dest="chain_length", default="100000000", help="Number of states to run the MCMC chain for. Default=100m")
    gen_group.add_argument("--log-every", dest="log_every", default="10000", help="How often to log tree and log files. Default=10,000")
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

    config = core_funcs.add_bools_to_config(config, args.multi_tree, args.fixed_tree, args.starting_tree, args.dta, args.glm, args.continuous_phylogeog)

    print(config)

    if args.fasta:
        config["taxa"], config["fasta"] = core_funcs.parse_fasta(args.fasta, args.codon_partitioning)
    else:
        config["taxa"] = core_funcs.get_taxa_no_fasta(args.id_file, args.id_file_dir, args.fixed_tree_file, config)
        config["fasta"] = False

    if config["fixed_tree"] or config["starting_tree"]:
        config["tree_name"], config["tree_file_dict"], config["newick_dict"] = core_funcs.fixed_tree_parsing(args.fixed_tree_file, args.starting_tree_file, args.fixed_tree_dir, config)
    else:
        config["tree_name"] = "tree1"
        config["tree_file_dict"] = False
        config["newick_dict"] = False

    if config["dta"]:
        config["ambiguities"] = trait_funcs.parse_ambiguities(args.ambiguities) #only for one trait atm
        config["traits"],config["trait_index"], config["all_trait_options"], config["trait_dict"], config["options_per_tree"] = trait_funcs.parse_discrete_traits(args.traits, args.discrete_trait_file,  args.trait_location_in_name, args.trait_delimiter, config)
    else:
        config["trait_index"] = False
        config["all_trait_options"] = False
        config["options_per_tree"] = False

    if config["glm"]: #needs dta bool to be true
        config["trait_to_predictor"], config["re_matrices"], config["bin_probs"] = glm_funcs.run_glm_functions(args.symmetric_predictor_dir, args.predictor_info_file, args.asymmetric_predictor_file, config)
    else:
        config["trait_to_predictor"] = False
        config["re_matrices"] = False
        config["bin_probs"] = False


    if config["continuous_phylogeog"]:
        #this template should be generalised to be fixed tree or not fixed, and multi/not multi
        config["traits"], config["trait_dict"], config["overall_trait"] = trait_funcs.continuous_phylogeography_processing(args.continuous_trait_file)
    else:
        config["overall_trait"] = False

    if not config["continuous_phylogeog"] and not config["dta"]:  
        config["traits"] = False
        config["trait_dict"] = False

    
    config["growth_model"] = args.growth_model
    if config["growth_model"] == "skygrid":
        config["gridpoints"] = int(args.sg_gridpoints)
        config["cutoff"] = args.sg_cutoff


    ##general options
    config["chain_length"] = args.chain_length
    config["log_every"] = args.log_every
    config["file_stem"] = args.file_stem
    config["template"] = args.template

    # for k,v in config.items():
    #     print(f'{k}: {v}')

    mytemplate = Template(filename=config["template"], strict_undefined=True)

    buf = StringIO()

    ctx = Context(buf, config=config)

    try:
        mytemplate.render_context(ctx)
    except:
        traceback = RichTraceback()
        for (filename, lineno, function, line) in traceback.traceback:
            print("File %s, line %s, in %s" % (filename, lineno, function))
            print(line, "\n")
        print("%s: %s" % (str(traceback.error.__class__.__name__), traceback.error))

    f = open(f"{config['file_stem']}.xml", 'w')
    # f.write(mytemplate.render(config=config))
    f.write(buf.getvalue())
    f.close()


if __name__ == '__main__':
    main()



