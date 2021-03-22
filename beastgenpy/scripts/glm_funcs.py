##Functions for making dictionaries containing predictors and other things required for GLMs

from Bio import SeqIO
from geopy import distance
from collections import defaultdict
import csv
import numpy as np
import statistics

def make_vector(matrix, dim):
    indices = np.triu_indices(dim,1) #get indices of top triangle, moved right by one to leave out diagonal
    upper = matrix[indices] #Create matrix of top triangle
    lowerprep = np.transpose(matrix) #Transpose matrix, so that lower triangle will be read in the desired order
    lower = lowerprep[indices] #Use same indices, but gets lower triangle of original matrix this time

    #turn arrays into strings so there are no commas 
    below_diag = np.array2string(lower, formatter = {'float':lambda lower: "%.1f" % lower}) 
    above_diag = np.array2string(upper, formatter = {'float': lambda upper: "%.1f" % upper})

    vector = above_diag + below_diag #make one vector that is the top triangle read horizontally, and lower triangle read vertically
    vector = vector.replace("]["," ").replace("[","").replace("]","")

    return vector

def process_info_file(info_file):

    standardised_transformed_list = []
    with open(info_file) as f:
        data = csv.DictReader(f)
        for l in data:
            if l["log_transformed_and_standardised"] == "True":
                standardised_transformed_list.append(l['predictor'])

    return standardised_transformed_list

def process_directional_predictors(trait_name, predictor_file, standardised_transformed_list):

    ##process folder, make sure 
    #format needs to be a column with a trait value and then one column per predictor value. Separate csvs for different traits
    predictor_list = []
    predictor_dict = defaultdict(dict)
    with open(predictor_file) as f:
        data = csv.DictReader(f)
        headers = data.fieldnames
        for i in headers:
            if i != trait_name:
                predictor_list.append(i)

        for predictor in predictor_list:
            inner_dict = {}
            for line in data:
                inner_dict[line[trait_name]] = float(line[predictor])

            if predictor in standardised_transformed_list:
                intermediate = {}
                for key,value in inner_dict.items():
                    intermediate[key] = np.log(values)
                new_inner = standardise(intermediate)
            else:
                new_inner = inner_dict
            
            predictor_dict[predictor] = new_inner


    sorted_predictor_values = defaultdict(list)
    for predictor, values_set in predictor_dict.items():
        ordered_values = sorted(values_set.items())
        for i in ordered_values:
            sorted_predictor_values[predictor].append(i[1])

    ## now going to make a matrix ##
    predictor_dict = {}
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

        predictor_dict[key_from] = from_value
        predictor_dict[key_to] = to_value

    return predictor_dict

def process_oneway_predictors(trait_name, trait_options, predictor_file):
           
    option_list = sorted(trait_options)
    dim = len(option_list)
    option_position = {}
    option_tups = []
    
    for i,option in enumerate(option_list):
        option_position[option] = i
        for option2 in option_list:
            option_tups.append((option, option2))

    tup_dict = {}
    for tup in option_tups:
        tup_dict[tup] = ""

    with open(predictor_file) as f:
        data = csv.DictReader(f)
        for l in data:
            first_ele = l[trait_name]
            for option in option_list:
                key = (first_ele, option)
                tup_dict[key] = l[option]
    
    matrix = np.zeros((dim, dim)) 
    for tup, value in tup_dict.items():
        location_1 = option_position[tup[0]]
        location_2 = option_position[tup[1]]

        matrix[location_1, location_2] = value

    vector = make_vector(matrix, dim)
    
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

def calculate_binomial_likelihoods(trait_options_dict):

    bin_probs  = {}
    
    for trait, options in trait_options_dict.items():
        n = len(options)
        p = 1 - (0.5**(1/n))
        bin_probs[trait] = p
    
    return bin_probs

def standardise(dictionary):
    
    standardised = {}

    sd = np.stdev(dictionary.values())
    mean = np.mean(dictionary.values())
    
    for key, value in dictionary.items():
        standardised[key] = ((value - mean)/sd)

    return standardised

def loop_for_processing(actual_predictor_dir, info_file, directional_file, trait_to_predictor):

    for f in os.listdir(actual_predictor_dir):

        if f != info_file and f != directional_file:

            predictor_name = f.strip(".csv")
            trait_to_predictor[trait_name][predictor_name] = process_oneway_predictors(trait_name, all_trait_options[trait_name], f)

