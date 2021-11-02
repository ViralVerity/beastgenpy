##Functions for making dictionaries containing predictors and other things required for GLMs

from Bio import SeqIO
from geopy import distance
from collections import defaultdict
import csv
import numpy as np
import statistics
import os

def run_glm_functions(predictor_dir_input, predictor_info_file, asymmetric_file, config):
#has to have predictors dir present if glm is true
#errors to do with the structure of the predictor dir

    re_matrices = make_twoway_REmatrices(config["all_trait_options"])
    bin_probs = calculate_binomial_likelihood(config["all_trait_options"])
    
    trait_to_predictor = defaultdict(dict)
    if len(config["traits"]) == 1:
        trait_to_predictor = glm_funcs.loop_for_processing(predictor_dir_input, predictor_info_file, asymmetric_file, config["traits"][0], trait_to_predictor, config["all_trait_options"])

    else:
        for trait_name in config["traits"]:
            predictor_dir = os.path.join(predictor_dir_input,trait)
            trait_to_predictor = glm_funcs.loop_for_processing(predictor_dir, predictor_info_file, asymmetric_file, trait_name, trait_to_predictor, config["all_trait_options"])

    return trait_to_predictor, re_matrices, bin_probs


def process_info_file(info_file):

    standardised_transformed_list = []
    with open(info_file) as f:
        data = csv.DictReader(f)
        for l in data:
            if l["log_transformed_and_standardised"].upper() == "TRUE":
                standardised_transformed_list.append(l['predictor'])

    return standardised_transformed_list

def process_asymmetric_predictors(trait_name, predictor_file, standardised_transformed_list):

    #format needs to be a column with the trait name (e.g. district in a phylogeography) and then one column per predictor. 
    #One csv for all the predictors of one trait.

    predictor_list = []
    predictor_dict = defaultdict(dict)

    with open(predictor_file) as f:
        data = csv.DictReader(f)
        headers = data.fieldnames
        for i in headers:
            if i != trait_name:
                predictor_list.append(i)
                
        for line in data:
            for predictor in predictor_list:
                predictor_dict[predictor][line[trait_name]] = float(line[predictor])
                
    final_predictor_dict = defaultdict(dict)
                
    for predictor, inner_dict in predictor_dict.items():
                
        if predictor in standardised_transformed_list:
            intermediate = {}
            for key,value in inner_dict.items():
                intermediate[key] = np.log(value)            
            new_inner = standardise(intermediate)
        else:
            new_inner = inner_dict

        final_predictor_dict[predictor] = new_inner

    sorted_predictor_values = defaultdict(list)
    for predictor, values_set in final_predictor_dict.items():
        ordered_values = sorted(values_set.items())
        for i in ordered_values:
            sorted_predictor_values[predictor].append(i[1])

    ## now going to make a matrix ##
    matrix_dict = {}
    for predictor, values in sorted_predictor_values.items():
        dim = len(values)
        frommatrix = np.zeros((dim,dim)) 
        tomatrix = np.zeros((dim,dim))

        key_from = f'{predictor}_origin'
        key_to = f'{predictor}_destination'

        for i, option in enumerate(values):
            frommatrix[i,:] = option
            tomatrix[:,i] = option
        
        from_value = make_vector(frommatrix, dim)  
        to_value = make_vector(tomatrix,dim)

        matrix_dict[key_from] = from_value
        matrix_dict[key_to] = to_value

    return matrix_dict

def process_symmetric_predictors(predictor_name, trait_name, trait_options, std_trans, predictor_file):
           
    option_list = sorted(trait_options)
    dim = len(option_list)
    option_position = {}
    option_tups = []

    for i,option in enumerate(option_list):
        option_position[option] = i
        for option2 in option_list:
            if option != option2:
                option_tups.append((option, option2))

    tup_dict = {}
    for tup in option_tups:
        tup_dict[tup] = ""

    with open(predictor_file) as f:
        data = csv.DictReader(f)
        headers = data.fieldnames
        for l in data:
            if trait_name in headers:
                first_ele = l[trait_name]
            elif predictor_name in headers:
                first_ele = l[predictor_name]
            else:
                first_ele = l[""]
            
            for option in option_list:
                key = (first_ele, option)
                if first_ele != option:
                    if std_trans:
                        tup_dict[key] = np.log(float(l[option]))
                    else:
                        tup_dict[key] = l[option]

    if std_trans:
        tup_dict = standardise(tup_dict)

    matrix = np.zeros((dim, dim)) 

    for tup, value in tup_dict.items():
        location_1 = option_position[tup[0]]
        location_2 = option_position[tup[1]]

        matrix[location_1, location_2] = value

    vector = make_vector(matrix, dim)

    return vector

def make_vector(matrix, dim):
    indices = np.triu_indices(dim,1) #get indices of top triangle, moved right by one to leave out diagonal
    upper = matrix[indices] #Create matrix of top triangle
    lowerprep = np.transpose(matrix) #Transpose matrix, so that lower triangle will be read in the desired order
    lower = lowerprep[indices] #Use same indices, but gets lower triangle of original matrix this time

    #turn arrays into strings so there are no commas 
    below_diag = np.array2string(lower, formatter = {'float':lambda lower: "%.1f" % lower}) 
    above_diag = np.array2string(upper, formatter = {'float': lambda upper: "%.1f" % upper})

    vector = above_diag + below_diag #make one vector that is the top triangle read horizontally, and lower triangle read vertically
    vector = vector.replace("]["," ").replace("[","").replace("]","").replace("\n","")

    return vector


def make_twoway_REmatrices(trait_options_dict):

    re_matrices = defaultdict(dict)

    for trait, options in trait_options_dict.items():
        trait_rand_design = {}

        #Create nxn matrices of zeros ready for populating
        dim = len(options)
        frommatrix = np.zeros((dim,dim)) 
        tomatrix = np.zeros((dim,dim))

        #These loops make from and to matrices for each trait, and then passes them through the function to get random effect vector
        for i, option in enumerate(options):
            key_from = f'{option}_from'
            key_to = f'{option}_to'
            if i == 0:
                frommatrix[i] = 1.0
                frommatrix[i,i] = 0.0

                tomatrix[:,i] = 1.0
                tomatrix[i,i] = 0.0
                
            else:
                frommatrix[i] = 1.0
                frommatrix[i-1] = 0.0
                frommatrix[i,i] = 0.0
                
                tomatrix[:,i] = 1.0
                tomatrix[:,i-1] = 0.0
                tomatrix[i,i] = 0.0
 
            trait_rand_design[key_from] = make_vector(frommatrix, dim)
            trait_rand_design[key_to] = make_vector(tomatrix,dim)

        re_matrices[trait] = trait_rand_design

    return re_matrices

def calculate_binomial_likelihood(trait_options_dict):

    bin_probs  = {}
    
    for trait, options in trait_options_dict.items():
        n = len(options)
        p = 1 - (0.5**(1/n))
        bin_probs[trait] = p
    
    return bin_probs

def standardise(dictionary):

    standardised = {}

    sd = np.std(list(dictionary.values()))
    mean = np.mean(list(dictionary.values()))
    
    for key, value in dictionary.items():
        standardised[key] = ((value - mean)/sd)

    return standardised

def loop_for_processing(actual_predictor_dir, info_file, asymmetric_file, trait_name, trait_to_predictor, all_trait_options):

    standardised_transformed_list = process_info_file(info_file)

    if actual_predictor_dir not in asymmetric_file:
        asymmetric_file = os.path.join(actual_predictor_dict, asymmetric_file)
    if actual_predictor_dir not in info_file:
        info_file = os.path.join(actual_predictor_dict, info_file)

    for f in os.listdir(actual_predictor_dir):
        pred_file = os.path.join(actual_predictor_dir, f)
        if pred_file != info_file and pred_file != asymmetric_file and pred_file.endswith(".csv"):
            if "/" in pred_file:
                predictor_name = pred_file.split("/")[-1].strip(".csv")
            else:
                predictor_name = pred_file.strip(".csv")
            if predictor_name in standardised_transformed_list:
                trait_to_predictor[trait_name][predictor_name] = process_symmetric_predictors(predictor_name, trait_name, all_trait_options[trait_name], True, pred_file)
            else:
                trait_to_predictor[trait_name][predictor_name] = process_symmetric_predictors(predictor_name, trait_name, all_trait_options[trait_name], False, pred_file)
        elif pred_file == asymmetric_file:
            trait_to_predictor[trait_name] = process_asymmetric_predictors(trait_name, asymmetric_file, standardised_transformed_list)

    return trait_to_predictor



###These functions haven't been integrated yet but will be useful when they have###
def pseudocount(inputfile, col_header):        
    pc_dict = {}
    with open(inputfile, "r") as f:
        r = csv.DictReader(f)
        data = [i for i in r]
        for loc in data:
            relevant_col = loc[col_header]
            pc_dict[loc["adm2"]] = relevant_col+1
            
    return pc_dict


def circle_distance(fw, loc1, loc2):
    d = distance.great_circle(centroid_dict[loc1], centroid_dict[loc2]).km
    fw.write(str(loc1) + "," + str(loc2) + "," + str(d) + '\n')

def get_centroids(map_file):

    locations = gp.read_file(map_file)

    locations['centroids'] = locations.centroid

    centroid_dict = {}
    for place, centroid in zip(locations["NAME_2"], locations["centroids"]):
        centroid_dict[place.upper().replace(","," ").replace(" ","_")] = ((centroid.x, centroid.y))

    location_pair_dict = defaultdict(list)

    for location in centroid_dict.keys():
        for location2 in centroid_dict.keys():
            location_pair_dict[location].append(location2)

    fw = open("distance_matrix.csv", 'w')
    for k,v in location_pair_dict.items():
        for item in v:
            circle_distance(fw, k, item)
    fw.close()

    return len(locations)

def write_ambiguity_codes(xml_file, ambig_file):

    with open(ambig_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split(",")
            adm2 = toks[0]
            ambigs = toks[1].replace("|"," ")
            xml_file.write(f'<ambiguity code="{adm2}" states="{ambigs}"/>\n')


def random_matrix_prolif(location, position):
    vector = np.zeros((dim,), dtype=int)
    vector[position] = 1
    line = '<parameter id="' + location + '_rand" value="' + str(vector) + '"/>' 
    line.replace("[", "")
