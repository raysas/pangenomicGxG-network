'''
This module contains the functions to construct the network from the data.  
Would need to start from the following data:  
- list of nodes (genes)
- GxS presence matrix
- phenotype matrix  
The goal is to compute the log odds of co-occurrence of genes in the presence matrix, and then use this to construct the network.

LOR of co-occurrence of genes i and j: log( a * d / b * c )  
where:
- a = number of R samples where both i and j are present
- b = number of R samples where they are not both present (1-a)
- c = number of S samples where both i and j are present
- d = number of S samples where they are not both present (1-c)

The network will be constructed as a weighted graph, where the weight of the edge between i and j is the LOR of co-occurrence of i and j.  
We will compute LOR for all possible pairs between the top genes, and take a threshold value to establish an edge between two genes.

Hence, the following attributes and functions are needed:  

Attributes:
-----------
- top_genes_file_path
- presence_path
- pheno_path

NOTE: add species and drug to the attributes so taht the user only sets these and everything will follow them

Functions:
----------
- Getters
- Setters

'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import os

os.chdir(os.path.expanduser('~/capstone-project'))


###################################### FUNCTIONS to get attributes ########################################################

def generate_phenotype_df(pheno_path:str, presence_df:pd.DataFrame=None)->pd.DataFrame:
    '''
    takes in a path to a phenotype file to make a df out of it
    if a presence_df is provided, it will filter the pheno_df to only include samples that are in the presence_df

    NOW TAHT TEH PHENO FILES ARE PROCESSED TO INCLUDE ONLY PRESENT GENOME IDS, THIS ONLY TRANSFORMS INTO DF

    param:
    ----
        - pheno_path: str
        - presence_df: pd.DataFrame (optional)

    return:
    ----
        - pheno_df: pd.DataFrame
    '''
    pheno_df = pd.read_csv(pheno_path, index_col=0)

    # -- already processed, this will give errors
    # if presence_df is not None:
    #     samples_presence= presence_df.columns
    #     mask = pheno_df.index.astype(str).isin(samples_presence)
    #     pheno_df = pheno_df.loc[mask]

    return pheno_df

def get_gene_list(file_path:str)->list:
    genes=[]
    if os.path.exists(file_path):
        with open(top_1000_nodes_path, 'r') as f:
            gene_list = f.readlines()
            genes = [line.strip() for line in gene_list]
    return genes


def generate_possible_gene_pairs(gene_list:list)->list:
    '''
    takes a list of genes and generates all possible gene pairs from the list

    param:
    ------
        - gene_list: list

    return:
    ------
        - gene_pairs: list
    '''
    gene_pairs=[]
    for i in range(len(gene_list)):
        for j in range(i+1, len(gene_list)):
            gene_pairs.append((gene_list[i], gene_list[j]))

    return gene_pairs

###################################### ATTRIBUTES ######################################################################

# -------------------- presence_path --------------------------------
presence_df=pd.DataFrame()
presence_path="data/temp/X_Escherichia_coli_GxS_filtered_hypo_unknowns_freq5.csv" #add the method that generates presence to this script
if os.path.exists(presence_path):
    presence_df= pd.read_csv(presence_path, index_col=0)
else:
    print("Please make sure you set the presence_path to the correct path")
    

# -------------------- pheno_path --------------------------------
pheno_df=pd.DataFrame()
pheno_path="metadata/Escherichia_coli/Escherichia_coli_trimethoprim.csv"
if os.path.exists(pheno_path):
    pheno_df= generate_phenotype_df(pheno_path, presence_df)
else:
    print("Please make sure you set the pheno_path to the correct path")

# -------------------- top_genes_file_path --------------------------------
top_1000_nodes_path = f'data/temp/SVM_Escherichia_coli_trimethoprim_top_1000_genes.txt'
gene_list = get_gene_list(top_1000_nodes_path)

# --------- one extra attribute to memoize the gene pairs ----------------
gene_pairs=generate_possible_gene_pairs(gene_list)


##############################################################################################################################

###################################### Getters & Setters ########################################################

def get_presence_path()->str:
    return presence_path

def get_presence_df()->pd.DataFrame:
    return presence_df

def set_presence_path(path:str):
    global presence_path
    presence_path=path
    global presence_df
    presence_df= pd.read_csv(presence_path, index_col=0)

def get_pheno_path()->str:
    return pheno_path

def get_pheno_df()->pd.DataFrame:
    return pheno_df

def set_pheno_path(path:str):
    global pheno_path
    pheno_path=path
    global pheno_df
    pheno_df= generate_phenotype_df(pheno_path, presence_df)

def get_top_genes_file_path()->str:
    return top_1000_nodes_path

def set_top_genes_file_path(path:str):
    global top_1000_nodes_path
    top_1000_nodes_path=path
    global gene_list
    gene_list=get_gene_list(top_1000_nodes_path)

def get_gene_pairs()->list:
    return gene_pairs

def set_gene_pairs(gene_list:list):
    global gene_pairs
    gene_pairs=generate_possible_gene_pairs(gene_list)

##############################################################################################################################

############################################### the good stuff ############################################################

def extract_genes_from_pairs(gene_pairs:list)->list:
    '''
    takes a list of gene pairs and returns a list of all genes in the gene pairs
    (opposite of the previous function)

    param:
    ------
        - gene_pairs: list of tuples

    return:
    ------
        - genes: list
    '''
    genes=[]
    for pair in gene_pairs:
        gene1=pair[0]
        gene2=pair[1]
        genes.append(gene1)
        genes.append(gene2)

    genes= list(set(genes))
    return genes

def split_matrix_by_phenotype(unlabeled_presence_df:pd.DataFrame, pheno_df: pd.DataFrame)->(pd.DataFrame, pd.DataFrame):
    '''
    takes a gene absence presence dataframe (n x m) with the samples classification and splits it to 2 dataframes:
        - one that represents all R samples (samples found in R_list)
        - one that represents all S samples (found in S_list)
    and returns them in this order: R then S

    param:
    ------
        - unlabeled_presence_df: (pd.DataFrame) the data frame output of data_utils.get_gene_presence_matrix
        - R_list: (list) samples genome_id that have a R (1) phenotype
        - S_list: (list) samples having the 0 phenotype
    
    return:
    ------
        - R_df: (pd.DataFrame) data frame where columns are only for R samples
        - S_df: (pd.DataFrame) df where cols are only for S samples

    '''


    # make sure pheno_df doesnthave more samples than those specified in the presence_df, subsample it otherwise
    samples_presence= unlabeled_presence_df.columns
    mask = pheno_df.index.astype(str).isin(samples_presence)
    pheno_df = pheno_df.loc[mask]

    # get the list of R and S samples:
    R=[];S=[];U=[]
    for sample in pheno_df.index:

        if pheno_df.loc[sample].values[0]==1:
            sample=str(sample)
            R.append(str(sample))
        elif pheno_df.loc[sample].values[0]==0:
            sample=str(sample)
            S.append(str(sample))
        else:
            sample=str(sample)
            U.append(str(sample))

    # All we care for care R and S which designate the list of resistant and susceptible samples respectively

    R_df = unlabeled_presence_df[R]
    S_df = unlabeled_presence_df[S]

    #make a subdf of all columns that are in S of the presence_df
    

    return R_df, S_df

def split_samples_list_by_phenotype(pheno_df:pd.DataFrame)->(list,list):
    '''
    takes a dataframe of the phenotypes and splits the samples into 2 lists of samples

    param:
    ------
        - pheno_df: pd.DataFrame
    return:
    ------
        - R_list: list
        - S_list: list
    '''
    # get the list of R and S samples:
    R=[];S=[];U=[]
    for sample in pheno_df.index:

        if pheno_df.loc[sample].values[0]==1:
            sample=str(sample)
            R.append(str(sample))
        elif pheno_df.loc[sample].values[0]==0:
            sample=str(sample)
            S.append(str(sample))
        else:
            sample=str(sample)
            U.append(str(sample))

    return R,S

# R,S=split_samples_list_by_phenotype(pheno_df)

def generate_RS_presence_counts(R:pd.DataFrame, S:pd.DataFrame)->pd.DataFrame:
    '''
    Takes the presence matrix of R samples and S samples and returns a dataframe of the count of genes present and absent in each group of samples.

    *This will memoize the entries needed for the log odds ratio computation.*

    
    e.g., output (w\o the log odds):

    | Gene       | R_present | R_absent | S_present | S_absent |
    |------------|-----------|----------|-----------|----------|
    | group_1001 | 96        | 0        | 187       | 0        |
    | tig        | 96        | 0        | 187       | 0        |
    | legF_1     | 96        | 0        | 187       | 0        |

    param:
    ----------------
        - R: (pd.DataFrame) presence matrix of R samples (output of get_subdf_cols)
        - S: (pd.DataFrame) presence matrix of S samples

    return:
    ----------------
        - new_df: (pd.DataFrame) dataframe of the count of genes present and absent in each group of samples.

    NOTE: it will perform the 0.5 correction - If for one gene any of the counts is 0, it will add 0.5 to all counts for that gene.

    '''
    R_present=R.sum(axis=1)
    R_absent=R.shape[1] - R_present
    S_present=S.sum(axis=1)
    S_absent=S.shape[1]-S_present


    new_df=pd.DataFrame({'R_present':R_present, 'R_absent':R_absent, 'S_present':S_present, 'S_absent':S_absent})

    #the 0.5 correction:
    for index in new_df.index:
        row=list(new_df.loc[index])
        if 0 in row: #check if any of the cols have value 0
            new_df.loc[index]=new_df.loc[index]+0.5

    return new_df

# RS_counts_df=generate_RS_presence_counts(R_df, S_df)

def compute_gene_log_odds_ratio(RS_counts_df:pd.DataFrame)->pd.DataFrame:
    '''
    takes the RS_counts_df and computes the log odds ratio for each gene in the dataframe

    param:
    ------
        - RS_counts_df: pd.DataFrame


    return:
    ------
        - log_odds_df: pd.DataFrame
    '''
    #get the total number of samples in each group
    n_R=len(R); n_S=len(S)


    R_present=RS_counts_df['R_present']
    R_absent=RS_counts_df['R_absent']
    S_present=RS_counts_df['S_present']
    S_absent=RS_counts_df['S_absent']

    log_odds=np.log((R_present/R_absent)/(S_present/S_absent))

    # print(log_odds)

    df=pd.DataFrame({"log_odds":log_odds})
    return df

def get_cooccurence_matrix(presence:pd.DataFrame)->pd.DataFrame:
    '''
    takes a binary GxS presence matrix and returns a GxG co-occurence matrix C 
    Cij= number of samples where gene i and gene j are both present
    (can be performed by matrix multiplication - trace on paper to understand why this works)

    input:
        - df: (pd.DataFrame) binary matrix of genes and samples
    output:
        - G: (nx.Graph) network of gene-gene presence
    '''
    gene_copresence= presence.dot(presence.T)
    return gene_copresence

def generate_RS_copresence_count(R_cooccurence:pd.DataFrame, S_cooccurence:pd.DataFrame, numR:int, numS:int, gene_pairs:list)->pd.DataFrame:
    '''
    Takes the co-occurence matrix of R samples and S samples and returns a dataframe of the count of genes copresent in each group of samples.
    The rows represent all gene-gene pairs and 4 columns are present: R_present, R_absent, S_present, S_absent

    gene_pairs is a list of tuples!

    *This will memoize the entries needed for the log odds ratio computation for gene-gene interaction*

    
    e.g., output (w\o the log odds):

    | Gene       | R_present | R_absent | S_present | S_absent |
    |------------|-----------|----------|-----------|----------|
    | G1-G2      | 96        | 0        | 187       | 0        |
    | G1-G3      | 96        | 0        | 187       | 0        |
    | G2-G3      | 96        | 0        | 187       | 0        |

    This example is without the 0.5 correction:  
        To ensure that the log odds ratio is defined, we will add 0.5 to all counts for each gene-gene pair if any of the counts is 0.

    param:
    ----------------
        - R_cooccurence: (pd.DataFrame) co-occurence matrix of R samples (output of get_cooccurence_matrix)
        - S_cooccurence: (pd.DataFrame) co-occurence matrix of S samples
        - numR: (int) number of R samples
        - numS: (int) number of S samples
        - gene_pairs: (list) list of gene-gene pairs that have a |correlation| of 0.6 or more (list of edges)

    return:
    ----------------
        - new_df: (pd.DataFrame) dataframe of the count of genes present in each group of samples.

    NOTE: it will perform the 0.5 correction - If for one gene any of the counts is 0, it will add 0.5 to all counts for that gene.
    '''


    df=pd.DataFrame(columns=['R_present', 'R_absent', 'S_present', 'S_absent'])

    for gene_pair in gene_pairs:
        gene1=gene_pair[0]
        gene2=gene_pair[1]
        # print(gene1, gene2)

        R_present=R_cooccurence.loc[gene1, gene2]
        R_absent=numR - R_present
        S_present=S_cooccurence.loc[gene1, gene2]
        S_absent=numS - S_present

        #the 0.5 correction:
        if R_present==0 or R_absent==0 or S_present==0 or S_absent==0:
            R_present+=0.5; R_absent+=0.5; S_present+=0.5; S_absent+=0.5

        df.loc[gene1+", "+gene2]=[R_present, R_absent, S_present, S_absent]

        # print(gene_pair, df.loc[gene1+"-"+gene2])

    # print(df)

    return df


def compute_association_log_odds(RS_cooccurence: pd.DataFrame)->pd.DataFrame:
    '''
    takes a RS_cooccurence G^2 x 4 matrix and computes the log odds of all gene pairs
    The idea is to compute the cooccurence log odds of resistance: 
        log ( a * d / c * b )  
        such that
        
         - a: is the number of strains that have both genes and are resistant
         - b: is the number of strains that do not have both genes and are resistant
         - c: is the of strains that have both genes and are susceptible
         - d: is the number of strains that do not have both genes and are susceptible

    param: 
    -------
    - RS_cooccurence: pd.Dataframe, G^2 x 4 matrix that contains the counts of copresence and coabsences of all gene pairs
    
    return:
    --------
    - log_odds_df: pd.DataFrame, the log odds calc of each pair of gene

    _all logs are computable given that RS_cooccurence have gone through the 0.5 correction_

    '''
    R_copresent=RS_cooccurence['R_present']
    R_coabsent=RS_cooccurence['R_absent']
    S_copresent=RS_cooccurence['S_present']
    S_coabsent=RS_cooccurence['S_absent']
    indices=RS_cooccurence.index
    # print(indices)

    log_odds=np.log((R_copresent*R_coabsent)/(S_copresent*S_coabsent))

    # print(log_odds)

    df=pd.DataFrame({"log_odds":log_odds}, index=indices)

    return df



def compute_cooccurence_LOR(gene_pairs:list=gene_pairs, pheno_df:pd.DataFrame=pheno_df, presence_df:pd.DataFrame=presence_df )->pd.DataFrame:
    '''
    takes a phenotype dataframe and a presence dataframe and computes the log odds ratio for each gene pair in the gene_pairs list

    param:
    ------
        - gene_pairs: list
        - pheno_df: pd.DataFrame
        - presence_df: pd.DataFrame

    return:
    ------
        - log_odds_df: pd.DataFrame
    '''
    # -- this step is important ro reduce run time and memory usage: filtering out genes from presence_df that we wont need to compute LOR
    genes=extract_genes_from_pairs(gene_pairs)
    presence_df=presence_df.loc[genes]

    R,S=split_samples_list_by_phenotype(pheno_df)
    R_df,S_df=split_matrix_by_phenotype(presence_df, pheno_df)
    R_cooccurence = get_cooccurence_matrix(R_df)
    S_cooccurence = get_cooccurence_matrix(S_df)
    numR=len(R); numS=len(S)

    RS_copresence_df=generate_RS_copresence_count(R_cooccurence, S_cooccurence, numR, numS, gene_pairs)
    log_odds_gene_pairs=compute_association_log_odds(RS_copresence_df)

    return log_odds_gene_pairs

def construct_network(log_odds_df:pd.DataFrame, threshold:float=0.5)->nx.Graph:
    '''
    takes a dataframe of the log odds of gene pairs and constructs a network from it

    param:
    ------
        - log_odds_df: pd.DataFrame
        - threshold: float

    return:
    ------
        - G: nx.Graph
    '''
    G=nx.Graph()

    for index in log_odds_df.index:
        log_odds=log_odds_df.loc[index].values[0]
        if log_odds>threshold:
            genes=index.split(", ")
            gene1=genes[0]
            gene2=genes[1]
            G.add_edge(gene1, gene2, weight=log_odds)

    return G

def construct_network_with_negative_weights(log_odds_df:pd.DataFrame, threshold:float=0.5)->nx.Graph:
    '''
    takes a dataframe of the log odds of gene pairs and constructs a network from it  

    the network constructed will be signed, in this case meaning that edges will occur between:
        - 2 highly R co-occurring genes
        - 2 highly S co-occurring genes
        - 1 highly R co-occurring gene and 1 highly S co-occurring gene

    This can be useful since in negative interactions are important interactions too

    param:
    ------
        - log_odds_df: pd.DataFrame
        - threshold: float

    return:
    ------
        - G: nx.Graph
    '''
    G=nx.Graph()

    for index in log_odds_df.index:
        log_odds=log_odds_df.loc[index].values[0]
        if log_odds>threshold or log_odds<-threshold:
            genes=index.split(", ")
            gene1=genes[0]
            gene2=genes[1]
            G.add_edge(gene1, gene2, weight=log_odds)

    return G

def workflow(pheno_path=pheno_path, presence_path=presence_path, gene_list_path=top_1000_nodes_path, threshold=0.5):
    '''
    This function is the main function that will run the entire workflow of constructing the network from the data

    param:
    ------
        - pheno_path: str
        - presence_path: str
        - gene_list_path: str
        - threshold: float

    return:
    ------
        - G: nx.Graph
    '''
    set_pheno_path(pheno_path)
    set_presence_path(presence_path)
    set_top_genes_file_path(gene_list_path)

    gene_pairs=generate_possible_gene_pairs(gene_list)
    log_odds_df=compute_cooccurence_LOR(gene_pairs, pheno_df, presence_df)
    G=construct_network(log_odds_df, threshold)

    return G