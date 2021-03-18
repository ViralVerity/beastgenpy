##Functions for making dictionaries containing predictors and other things required for GLMs

from Bio import SeqIO
from geopy import distance
from collections import defaultdict
import csv
import numpy as np
import statistics

def REmatrices(matrix):
    indices = np.triu_indices(dim,1) #get indices of top triangle, moved right by one to leave out diagonal
    upper = matrix[indices] #Create matrix of top triangle
    lowerprep = np.transpose(matrix) #Transpose matrix, so that lower triangle will be read in the desired order
    lower = lowerprep[indices] #Use same indices, but gets lower triangle of original matrix this time

    #turn arrays into strings so there are no commas 
    below_diag = np.array2string(lower, formatter = {'float':lambda lower: "%.1f" % lower}) 
    above_diag = np.array2string(upper, formatter = {'float': lambda upper: "%.1f" % upper})

    vector = above_diag + below_diag #make one vector that is the top triangle read horizontally, and lower triangle read vertically

    return vector

def make_twoway_REmatrices(trait_dict):

    re_matrices = defaultdict(dict)

    for trait, options in trait_dict.itesm():
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
                from_value = REmatrices(frommatrix)

                tomatrix[:,i] = 1.0
                tomatrix[i,i] = 0.0
                to_value = REmatrices(tomatrix)
            else:
                frommatrix[i] = 1.0
                frommatrix[i-1] = 0.0
                frommatrix[i,i] = 0.0
                from_value = REmatrices(frommatrix)

                tomatrix[:,i] = 1.0
                tomatrix[:,i-1] = 0.0
                tomatrix[i,i] = 0.0
                to_value = REmatrices(tomatrix)

            trait_rand_design[key_from] = from_value
            trait_rand_design[key_to] = to_value

        re_matrices[trait] = trait_rand_design

    return re_matrices
