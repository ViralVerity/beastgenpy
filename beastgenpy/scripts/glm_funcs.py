from Bio import SeqIO
from geopy import distance
from collections import defaultdict
import csv
import numpy as np
import statistics


def pseudocount(inputfile, col_header):        
    pc_dict = {}
    with open(inputfile, "r") as f:
        r = csv.DictReader(f)
        data = [i for i in r]
        for loc in data:
            relevant_col = loc[col_header]
            pc_dict[loc["adm2"]] = relevant_col+1
            
    return pc_dict


def transform_data(dictionary):
    
    transformed = {}

    sd = np.stdev(dictionary.values())
    mean = np.mean(dictionary.values())
    
    for key, value in dictionary.items():
        transformed[key] = ((value - mean)/sd)

    return transformed


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

def REmatrices(matrix):
    indices = np.triu_indices(dim,1) #get indices of top triangle, moved right by one to leave out diagonal
    upper = matrix[indices] #Create matrix of top triangle
    lowerprep = np.transpose(matrix) #Transpose matrix, so that lower triangle will be read in the desired order
    lower = lowerprep[indices] #Use same indices, but gets lower triangle of original matrix this time

    #turn arrays into strings so there are no commas 
    below_diag = np.array2string(lower, formatter = {'float':lambda lower: "%.1f" % lower}) 
    above_diag = np.array2string(upper, formatter = {'float': lambda upper: "%.1f" % upper})

    vector = above_diag + below_diag #make one vector that is the top triangle read horizontally, and lower triangle read vertically

def make_REmatrices(dim):

    #Create nxn matrices of zeros ready for populating
    frommatrix = np.zeros((dim,dim)) 
    tomatrix = np.zeros((dim,dim))

    #These loops make from and to matrices for each location, and then passes them through the function to get random effect vector
    for i in range(dim):
        if i == 0:
            frommatrix[i] = 1.0
            frommatrix[i,i] = 0.0
            REmatrices(frommatrix)
        else:
            frommatrix[i] = 1.0
            frommatrix[i-1] = 0.0
            frommatrix[i,i] = 0.0
            REmatrices(frommatrix)


    for i in range(dim):
        if i == 0:
            tomatrix[:,i] = 1.0
            tomatrix[i,i] = 0.0
            REmatrices(tomatrix)
        else:
            tomatrix[:,i] = 1.0
            tomatrix[:,i-1] = 0.0
            tomatrix[i,i] = 0.0
            REmatrices(tomatrix)

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

    




