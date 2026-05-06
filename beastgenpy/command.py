import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import argparse
import os
import csv
from mako.template import Template
from mako.runtime import Context
from mako.exceptions import RichTraceback
from mako.template import Template
from mako.lookup import TemplateLookup
from io import StringIO
from collections import defaultdict
import datetime as dt

import glm_funcs as glm_funcs
import core_funcs as core_funcs
import trait_analysis_funcs as trait_funcs
import error_checks as error_checks

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False, description='make XML')

    tax_group = parser.add_argument_group('Taxa options')
    tax_group.add_argument("--alignment", help="comma separated list of fasta files containing alignment to analyse")
    tax_group.add_argument("--codon-partitioning", "-cp", dest="codon_partitioning", help="comma separated list of 1s and 0s for which alignments have codon partitioning")

    tree_group = parser.add_argument_group("Tree options")
    tree_group.add_argument("--empirical-tree", action="store_true", dest="empirical", help="perform analysis on empirical trees")
    tree_group.add_argument("--fixed-tree", action="store_true", dest="fixed_tree", help="Perform analysis on fixed tree")
    tree_group.add_argument("--starting-tree", action="store_true", dest="starting_tree", help="flag for adding a starting tree")
    tree_group.add_argument("--tree-file", dest="fixed_tree_file", help="File with single fixed/empirical/starting tree") 
    tree_group.add_argument("--tree-dir", dest='fixed_tree_dir', help="Directory containing multiple fixed/empirical/starting trees in newick format")

    population_model_group = parser.add_argument_group("population model tree priors")
    population_model_group.add_argument("--population-model", dest="population_model", default="skygrid")
    population_model_group.add_argument("--sg-cutoff", dest="sg_cutoff")
    population_model_group.add_argument("--sg-gridpoints", dest="sg_gridpoints")

    clock_model_group = parser.add_argument_group("clock models")
    clock_model_group.add_argument("--clock-model", dest="clock_model", default="relaxed")

    subs_model_group = parser.add_argument_group("subs models")
    subs_model_group.add_argument("--subs-model", dest="subs_model", default="gtr")

    trait_group = parser.add_argument_group("trait_analysis_group")
    trait_group.add_argument("--phylogeography", help="trait analysis, options are discrete or continuous")
    trait_group.add_argument("--trait-file", dest="trait_file", help="file containing values for discrete or continuous traits")
    trait_group.add_argument("--polygon-dir", help="directory with polygons for uncertainty estimation", dest="polygon_dir")

    trait_group.add_argument("--traits", help="Comma separated list of traits for discrete trait analysis")
    trait_group.add_argument("--trait-location-in-name", dest="trait_location_in_name", help="Information to pull trait value from taxon label. Should be provided as trait=position_in_name with different traits separated by commas. e.g.location=5, counting from 1")
    trait_group.add_argument("--trait-delimiter", dest="trait_delimiter", default="|", help="Separater in taxon label to pull trait value from taxon label. Default='|'")
    trait_group.add_argument("--ambiguities", help="tsv containing ambiguities for dta with two columns: one with header 'ambiguity' and one with 'options' containing a comma separated list")

    trait_group.add_argument("--glm", action="store_true", help="Flag to run a generalised linear model on discrete trait predictors")
    trait_group.add_argument("--predictor-info-file", dest="predictor_info_file", help="File containing which predictors should be logged and standardised in csv format. Headers must be 'predictor' and 'logged_transformed_and_standardised'") 
    trait_group.add_argument("--asymmetric-predictor-file", dest="asymmetric_predictor_file", help="File containing asymmetric i.e. one-way predictors for a trait in csv format. Headers must be the name of the trait, and then the predictor names") 
    trait_group.add_argument("--symmetric-predictor-dir", dest="symmetric_predictor_dir", help="Directory containing matrices for symmetrical predictors. Each file should have the first row and column as trait values, with the predictor value in the appropriate row/column combination.")
    trait_group.add_argument("--epoch", action="store_true", help="flag to make it an epoch analysis")
    trait_group.add_argument("--transition-times", dest="transition_times", help="comma separated list of values for epoch transition times in terms of years from most recent tip")
    
    
    gen_group = parser.add_argument_group('General options')
    gen_group.add_argument("--chain-length", dest="chain_length", default="100000000", help="Number of states to run the MCMC chain for. Default=100m")
    gen_group.add_argument("--log-every", dest="log_every", default="10000", help="How often to log tree and log files. Default=10,000")
    gen_group.add_argument("--file-stem", dest="file_stem", help="File stem for analysis")
    gen_group.add_argument("--verbose", help="print more info")
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

    config = core_funcs.add_bools_to_config(config, args.fixed_tree, args.empirical, args.starting_tree, args.glm, args.epoch, args.verbose)

    #pull information about sequences
    if args.alignment:
        config["sequence_info"] = core_funcs.parse_fasta(args.alignment, args.codon_partitioning)
    if args.fixed_tree:
        config = core_funcs.parse_fixed_trees(args.tree_file, args.tree_dir)
    if not config["sequence_info"]:
        sys.stderr.write("Need to provide either an alignment or fixed trees file/directory\n")
        sys.exit(-1)

    error_checks.check_dates_in_names(config)

    config["seq_to_tree"] = {}
    for name, info in config["sequence_info"].items():
        for seq_name in info["taxon_list"]:
            config["seq_to_tree"][seq_name] = name

    if args.starting_tree:
        config["sequence_info"] = core_funcs.parse_starting_trees(args.tree_file, args.tree_dir)

    if args.phylogeography:
        config = error_checks.check_phylogeog_value(config, args.phylogeography)

        if config["phylogeography"] == "discrete":
            config["ambiguities"] = trait_funcs.parse_ambiguities(args.ambiguities) #only for one trait atm
            config["traits"],config["trait_index"], config["all_trait_options"], config["trait_dict"], config["options_per_tree"] = trait_funcs.parse_discrete_traits(args.traits, args.trait_file,  args.trait_location_in_name, args.trait_delimiter, config)
        
            if config["glm"]: 
                config["trait_to_predictor"], config["re_matrices"], config["bin_probs"] = glm_funcs.run_glm_functions(args.symmetric_predictor_dir, args.predictor_info_file, args.asymmetric_predictor_file, config)
                config = glm_funcs.get_markov_counts(config)
            else:
                config["trait_to_predictor"] = False
                config["re_matrices"] = False
                config["bin_probs"] = False
        else:
            config["trait_index"] = False
            config["all_trait_options"] = False
            config["options_per_tree"] = False

        
        if config["phylogeography"] == "continuous":
            config["traits"], config["trait_dict"], config["overall_trait"] = trait_funcs.continuous_phylogeography_processing(args.trait_file)
            error_checks.check_seqs_present(config)
            if args.polygon_dir:
                error_checks.check_file(args.polygon_dir)
                config["uncertain_polygons"] = trait_funcs.sort_uncertain_polygons(args.polygon_dir)
                config["polygon_dir"] = args.polygon_dir
            else:
                config["uncertain_polygons"] = []
        else:
            config["overall_trait"] = False


    else:
        config["traits"] = False
        config["trait_dict"] = False
        config["phylogeography"] = False

    if config["epoch"]:
        config["transition_times"] = args.transition_times.split(",")
    else:
        config["transition_times"] = []

    
    config["population_model"] = args.population_model
    if config["population_model"] == "skygrid":
        config = error_checks.check_gp_cutoff(config, args.sg_cutoff, args.sg_gridpoints)
        if not config["gridpoints"]:
            config["gridpoints"] = int(args.sg_gridpoints)

    config["clock_model"] = args.clock_model
    config["subs_model"] = args.subs_model

    if args.file_stem:
        config["file_stem"] = args.file_stem
    else:
        config["file_stem"] = "beast"
    
    ##general options
    config["chain_length"] = args.chain_length
    config["log_every"] = args.log_every

    path_to_templates = os.path.join(thisdir, "templates")
    mylookup = TemplateLookup(directories=[path_to_templates])
    mytemplate = Template(filename=os.path.join(path_to_templates, "master_template.template"), uri="master_template.template", strict_undefined=True, lookup=mylookup)

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



