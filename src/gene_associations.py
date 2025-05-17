'''       
This module is to compute different statistical measures to associate genes with the phenotype.  
important considerations: setting a species and a drug of interest.

The main statistical associations will be:  
- Mutual Information
- Chi-Square test
- ANOVA test

The main ML models to infer feature importance will be:
- SVM ensemble

Attributes
-----------
- species: the species name
- drug: the drug of interest

Derived attributes:
--------------------
- presence_path: the path to the presence matrix
- pheno_path: the path to the phenotype file

- presence_matrix: the presence matrix of the species
- labeled_matrix: the labeled matrix of the species (maybe create a function for that)


'''
import os
import sys
import pandas as pd
import numpy as np

os.chdir(os.path.expanduser('~/capstone-project'))
sys.path.append('src')

import cluster_analysis

# ---------------------- validation functions (retrival and rank) ----------------------

def parse_ARGs_from_df(ARG_products, top_genes)-> dict:
    '''
    Takes in a list of ARG products and a list of top genes and returns a dict of ARGs that are in the top genes (with their rank)
    
    param:
    ------
    - ARG_products: ARG products
    - top_genes: of top genes

    return:
    -------
    - ARGs: dict, dict of ARGs that are in the top genes (key rank, value gene - maybe consider switching)
    '''
    ARGs = {}
    for gene in top_genes:
        for product in ARG_products:
            if product in gene:
                ARGs[top_genes.index(gene)+1]=gene # +1 to start ranking from 1
    return ARGs 

def get_ranked_ARGs_from_association(ass_score_dict:dict, products_list, n=100, sort_reverse=True):
    '''
    takes a dictionary of association score (key cluster and value the score) & the list of products we want to find their rank in a sorted list of te top 100 associations
    gets us a dict of key rank and value the ARG product; 
    sort_reverse is set to True when the higher value means higher rank, and to false when lower value means higher rank (the case with p values as values)

    param:
    ------
    - ass_score_dict: dict
    - products_list: list
    - n: int, 100 by default

    return:
    ---------
    - ranked_dict: dict
    '''
    top_100 = sorted(ass_score_dict, key=ass_score_dict.get, reverse=sort_reverse)[:n]
    for i in  range(len(top_100)):
        top_100[i]=cluster_analysis.transform_cluster_to_product(top_100[i])

    ranked_dict = parse_ARGs_from_df(products_list, top_100)
    return ranked_dict

# ---------------------- statistical tests and associations ----------------------
