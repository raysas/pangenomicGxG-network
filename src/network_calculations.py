'''
This script is intended to compute correlations, associations and other metrics for network features:

helper methods:
    - _get_subdf_cols: (pd.DataFrame, list) -> pd.DataFrame
        returns a sub dataframe with only the columns in the list
        used in _split_matrix_by_phenotype
    - _split_matrix_by_phenotype: (pd.DataFrame, list, list) -> pd.DataFrame, pd.DataFrame
        returns 2 dataframes, one for R samples and one for S samples
        used in log_odds
    - _label_samples: (pd.DataFrame, list) -> list, list, list
        returns 3 lists of resistant, susceptible and unclassfied samples
        used in log_odds
    - _create_RS_presence_stats_matrix: (pd.DataFrame, pd.DataFrame) -> pd.DataFrame
        returns a dataframe of the count of genes present and absent in each group of samples
        used in log_odds
    
main method:
    - log_odds: (str, str, str) -> pd.DataFrame
        returns the log odds ratio matrix of the genes
        used in main.py
'''

import numpy as np
import pandas as pd
import networkx as nx
from data_utils import *
import matplotlib.pyplot as plt

_get_subdf_cols= lambda df, col_list: df[col_list]

def _split_matrix_by_phenotype(unlabeled_presence_df:pd.DataFrame, R_list: list, S_list: list):
    '''
    takes a gene absence presence dataframe (n x m) with the samples classification and splits it to 2 dataframes:
        - one that represents all R samples (samples found in R_list)
        - one that represents all S samples (found in S_list)
    and returns them in this order: R then S

    param:
        - unlabeled_presence_df: (pd.DataFrame) the data frame output of data_utils.get_gene_presence_matrix
        - R_list: (list) samples genome_id that have a R (1) phenotype
        - S_list: (list) samples having the 0 phenotype
    
    return:
        - R_df: (pd.DataFrame) data frame where columns are only for R samples
        - S_df: (pd.DataFrame) df where cols are only for S samples

    _this is a helper function to create a contigency table for each gene_
    '''
    R_df=_get_subdf_cols(unlabeled_presence_df, R_list)
    S_df=_get_subdf_cols(unlabeled_presence_df, S_list)

    return R_df, S_df

def _label_samples(pheno_df: pd.DataFrame) :
    '''
    Takes a dataframe of phenotypes out of the csv file of a species-drug, a list of all samples and returns 3 lists of resistant, susceptible and unclassfied samples
    In the csv if 
    param:
        - pheno_df: (pd.DataFrame) dataframe of phenotypes
        - samples_list: (list) list of all samples  

    output:
        - R: (list) list of resistant samples
        - S: (list) list of susceptible samples
        - U: (list) list of unclassified samples
    '''
    R=[];S=[];U=[]
    for sample in pheno_df.index:

        if pheno_df.loc[sample].values[0]==1:
            sample=str(sample)
            if len(sample)==7:
                sample+='0'
            elif len(sample)==6:
                sample+='00'
            elif len(sample)==5:
                sample+='000'
            elif len(sample)==4:
                sample+='0000'
            R.append(str(sample))
        elif pheno_df.loc[sample].values[0]==0:
            sample=str(sample)
            sample=str(sample)
            if len(sample)==7:
                sample+='0'
            elif len(sample)==6:
                sample+='00'
            elif len(sample)==5:
                sample+='000'
            elif len(sample)==4:
                sample+='0000'
            S.append(str(sample))
        else:
            sample=str(sample)
            if len(sample)==7:
                sample+='0'
            elif len(sample)==6:
                sample+='00'
            elif len(sample)==5:
                sample+='000'
            elif len(sample)==4:
                sample+='0000'
            U.append(str(sample))


def _create_RS_presence_stats_matrix(R: pd.DataFrame, S: pd.DataFrame) -> pd.DataFrame:
    '''
    takes the presence matrix of R samples and S samples and returns a dataframe of the count of genes present and absent in each group of samples.

    *This will memoize the entries needed for the log odds ratio computation.*

    
    _e.g., output (w\o the log odds):_

    | Gene       | R_present | R_absent | S_present | S_absent |
    |------------|-----------|----------|-----------|----------|
    | group_1001 | 96        | 0        | 187       | 0        |
    | tig        | 96        | 0        | 187       | 0        |
    | legF_1     | 96        | 0        | 187       | 0        |

    param:
        - R: (pd.DataFrame) presence matrix of R samples (output of get_subdf_cols)
        - S: (pd.DataFrame) presence matrix of S samples

    return:
        - new_df: (pd.DataFrame) dataframe of the count of genes present and absent in each group of samples.


    '''
    R_present=R.sum(axis=1)
    R_absent=R.shape[1] - R_present
    S_present=S.sum(axis=1)
    S_absent=S.shape[1]-S_present

    #for each row that has all e

    log_odds=np.log((R_present/R_absent)/(S_present/S_absent))

    new_df=pd.DataFrame({'R_present':R_present, 'R_absent':R_absent, 'S_present':S_present, 'S_absent':S_absent, 'log_odds':log_odds})
    return new_df

def log_odds(pheno_path:str, Rtab_presence_matrix_path:str, remove_hypothetical:bool=False):
    '''
    Takes the path to the phenotypes csv, the species ids file and the gene presence/absence matrix and returns the log odds ratio matrix of the genes.

    param:
    -------
        - pheno_path: (str) path to the phenotypes csv
        - Rtab_presence_matrix_path: (str) path to the gene presence/absence matrix
        - remove_hypothetical: (bool) whether to remove hypothetical genes or not
        
    return:
    -------
        - log_odds_dict: (dict) dictionary of the log odds ratio of the genes

    The output can be further used as node attributes in the network.
    '''
    pheno=get_pheno_df(pheno_path)
    presence=get_gene_presence_matrix(Rtab_presence_matrix_path, remove_hypothetical)

    R,S,U=_label_samples(pheno)
    R_df, S_df=_split_matrix_by_phenotype(presence, R, S)

    stats_table=_create_RS_presence_stats_matrix(R_df, S_df)
    log_odds_table=stats_table['log_odds']

    log_odds_dict = log_odds_table.to_dict()

    return log_odds_dict


def get_network_stats(G:nx.Graph, network_name:str='test', weighted:bool=False)->pd.DataFrame:
    '''
    Takes a graph and returns the stats in a dataframe

    param:
        - G: (nx.Graph) the graph
        - weighted: (bool) whether the graph is weighted or not
        - network_name: (str) the name of the network

    return:
        - stats_df: (pd.DataFrame) dataframe of the stats
    '''
    components=nx.connected_components(G)
    largest_cc = max(components, key=len)
    largest_cc_subgraph = G.subgraph(largest_cc)
    number_of_components=len(components)
    communities=nx.community.louvain_communities(G)


    if weighted:
        stats = {
            "Network": network_name,
            "|N|": G.number_of_nodes(),
            "|V|": G.number_of_edges(),
            "Density": nx.density(G),
            "k": np.mean(list(dict(G.degree()).values())),
            "k weighted": np.mean(list(dict(G.degree(weight='weight')).values())),
            "cc": nx.average_clustering(G),
            "spath": nx.average_shortest_path_length(largest_cc_subgraph),
            "d": nx.diameter(largest_cc_subgraph),
            "Q": nx.community.modularity(G, communities),
            '# of components': number_of_components,
            "# of communities": len(communities),
        }
    else:
        stats = {
            "Network": network_name,
            "|N|": G.number_of_nodes(),
            "|V|": G.number_of_edges(),
            "Density": nx.density(G),
            "<k>": np.mean(list(dict(G.degree()).values())),
            "<cc>": nx.average_clustering(G),
            "<spath>": nx.average_shortest_path_length(largest_cc_subgraph),
            "d": nx.diameter(largest_cc_subgraph),
            "Q": nx.community.modularity(G,communities ),
            '# of components': number_of_components,
            "# of communities": len(communities)
        }
    
    stats_df = pd.DataFrame(stats, index=[0])
    stats_df.set_index('Network', inplace=True)

    return stats_df

def plot_degree_distribution(G:nx.Graph, title:str='', weighted:bool=False):
    '''
    Takes a graph and plots the degree distribution

    param:
        - G: (nx.Graph) the graph
        - weighted: (bool) whether we want the weighted deg dist
        - title: (str) info on the network

    return:
        - None
    '''
    if not weighted:
        degrees = [G.degree(n) for n in G.nodes()]
        plt.hist(degrees, bins=30,color='lightblue')
        plt.title(network_name)
        plt.show()
    else:
        degrees = [G.degree(n, weight='weight') for n in G.nodes()]
        plt.hist(degrees, bins=30, label=network_name)
        plt.title(network_name)
        plt.show()

def set_log_odds_node_attributes(G:nx.Graph, log_odds:pd.DataFrame):
    '''
    Takes a graph and sets the log odds as node attributes

    param:
        - G: (nx.Graph) the graph
        - log_odds: (pd.DataFrame) the log odds dataframe

    return:
        - None
    '''
    for gene in G.nodes():
        if gene in log_odds.index:
            G.nodes[gene]['log_odds'] = log_odds.loc[gene].values[0]
        else:
            G.nodes[gene]['log_odds'] = 0