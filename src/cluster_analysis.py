'''
Cluster Analysis Module
-----------------------
Responsible for extracting the cluster attributes from the pipeline output and setting them as node attributes in the networkx graph

The pipeline output consists of 4 files that we need:

    - *.fatsa.clstr file (cd hit output) that contains the cluster number and the sequences in all of the genomes
    - *_cluster_frequencies.csv file that contains the cluster number and the frequency of the cluster in the genomes
    - *_pangenome.csv file that contains the cluster number and the gene class of the cluster
    - *.fasta file (cd hit output) that contains the cluster number and the representative sequence of the cluster4

Attributes:
----------
    _output_path: one global (private) attribute that is the path to the pipeline output directory - changed using a setter `set_output_path(str)`

Functions: 
----------
 
note: even though there are many helper functions, kept public in case needed (so far)

* autoincrement_col(df:pd.DataFrame, col:str)->pd.DataFrame
* generate_files_paths(pipeline_output_path:str)->str, str, str, str
* get_cluster_frequency(freq_path:str)->pd.DataFrame
* get_cluster_representatives(clstr_path:str, save_csv:bool=False)->pd.DataFrame
* get_representative_products(clstr_fasta_path:str)->pd.DataFrame
* combine_cluster_product(clstr_rep_df:pd.DataFrame, rep_prod_df:pd.DataFrame)->pd.DataFrame
* get_cluster_pan_gene_class(pan_annot_path:str)->pd.DataFrame
* get_cluster_representatives(clstr_path:str, save_csv:bool=False)->pd.DataFrame
* get_cluster_attributes(clstr_path, freq_path, pan_annot_path, clstr_fasta_path)
* create_dicts(df:pd.DataFrame)->dict
* set_node_attributes_by_cluster(df:pd.DataFrame, G:nx.graph)->nx.Graph
* set_output_path(path:str)->str
* set_cluster_attributes(G:nx.Graph, pipeline_output=_output_path)->nx.Graph
* split_matrix_by_phenotype(unlabeled_presence_df:pd.DataFrame, pheno_df: pd.DataFrame)->(pd.DataFrame, pd.DataFrame)
* generate_RS_presence_counts(R:pd.DataFrame, S:pd.DataFrame)->pd.DataFrame
* compute_log_odds_ratio(RS_counts_df:pd.DataFrame)->pd.DataFrame
* get_cluster_resistance_LOR(presence_df:pd.DataFrame, pheno_df:pd.DataFrame)-> pd.DataFrame
* transform_cluster_to_product(cluster_name, clstr_product_df)
    
'''


import requests
import re
import os
import sys
import networkx as nx
import pandas as pd
import numpy as np

os.chdir(os.path.expanduser('~/capstone-project'))
sys.path.append("src")

# -------------- Attributes ----------------

species='Escherichia_coli'
_output_path=f'data/pangenome_pipeline_output/{species}'

# -- getters and setters

def set_output_path(path:str):
    '''
    sets the output path for the pipeline

    param:
    ------
    - path: str, path to the output directory
    '''
    global _output_path
    _output_path=path

    global clstr_path, freq_path, pan_annot_path, clstr_fasta_path
    clstr_path, freq_path, pan_annot_path, clstr_fasta_path= generate_files_paths(_output_path)


    return _output_path

def get_output_path():
    '''
    returns the output path for the pipeline
    '''
    return _output_path

def set_species(sp:str):
    '''
    sets the species for the pipeline, changing it will change the output path and the 4 files paths

    param:
    ------
    - species: str, species name
    '''
    global species
    species=sp
    global _output_path
    _output_path=f'data/pangenome_pipeline_output/{species}'

    global clstr_path, freq_path, pan_annot_path, clstr_fasta_path
    clstr_path, freq_path, pan_annot_path, clstr_fasta_path= generate_files_paths(_output_path)

    return species

def get_species():
    '''
    returns the species for the pipeline
    '''
    return species

# -------------- Helper Functions ----------------

def autoincrement_col(df:pd.DataFrame, col:str)->pd.DataFrame:
    '''
    Autoincrement a column by adding _# to it everytime its duplicated

    e.g., one -> one_1, one -> one_2, two -> two_1

    param:
    ------
    - df: pd.DataFrame, dataframe
    - col: str, column name to autoincrement

    return:
    -------
    - df: pd.DataFrame, dataframe with the column autoincremented

    p.s. the resultant dataframe will be have unqiue values in the unqiue - can be firther set as index if desired
    '''
    df[col] = df[col].add('_').add(df.groupby(col).cumcount().add(1).astype(str))
    return df


def generate_files_paths(pipeline_output_path):
    '''
    takes the output dir path of the pipeline for a species, and returns the 4 file paths inside of it needed to extract the data

    Parameters:
    -----------
    pipeline_output_path: str, path to the pipeline output directory

    Returns:
    --------
    clstr_path: str, path to the cluster file, cdhit output
    freq_path: str, path to the cluster frequency file
    pan_annot_path: str, path to the pangenome annotation file
    clstr_fasta_path: str, path to the cluster fasta file, aslo cd hit output
    '''

    pipeline_output_path=pipeline_output_path.rstrip('/')
    basename=os.path.basename(pipeline_output_path)

    clstr_path=pipeline_output_path+'/'+basename+'.fasta.clstr'
    freq_path=pipeline_output_path+'/'+basename+'_cluster_frequencies.csv'
    pan_annot_path=pipeline_output_path+'/'+basename+'_pangenome.csv'
    clstr_fasta_path=pipeline_output_path+'/'+basename+'.fasta'

    return clstr_path, freq_path, pan_annot_path, clstr_fasta_path

# -------------- deeclare new attributes ---------------
clstr_path, freq_path, pan_annot_path, clstr_fasta_path= generate_files_paths(_output_path)

# -- getters and setters

def get_clstr_path():
    return clstr_path

def get_freq_path():
    return freq_path

def get_pan_annot_path():
    return pan_annot_path

def get_clstr_fasta_path():
    return clstr_fasta_path

def set_clstr_path(path:str):
    global clstr_path
    clstr_path = path

def set_freq_path(path:str):
    global freq_path
    freq_path = path

def set_pan_annot_path(path:str):
    global pan_annot_path
    pan_annot_path = path

def set_clstr_fasta_path(path:str):
    global clstr_fasta_path
    clstr_fasta_path = path

#  ------------------------------------------------------------------------------------
    

def get_cluster_frequency(freq_path):
    '''
    Takes the path to the cluster frequency file and returns a dataframe with the cluster frequencies

    Parameters:
    freq_path: str

    Returns:
    freq_df: pandas.DataFrame
    '''

    freq_df=pd.read_csv(freq_path, index_col=0)
    freq_df.columns=['Cluster','frequency']
    freq_df.index=freq_df['Cluster']
    freq_df.drop('Cluster', axis=1, inplace=True)

    return freq_df

def get_cluster_representatives(clstr_path:str, save_csv:bool=False)->pd.DataFrame:
    '''
    Take a CD-HIT clstr output file and return a dictionary with cluster number as key and representative sequence as value.
    Optionally saves the dictionary as a csv file.

    parameters:
    ----------
    clstr_path: str, path to CD-HIT clstr file
    save_csv: bool, if True saves the dictionary as a csv file (optional)

    return:
    -------
    df: pd.DataFrame, data frame with cluster number as key and representative sequence as value

    description:
    ------------

    CD-HIT cluster file parser:
        .clstr file format is:

        >Cluster 489
        0	552aa, >fig|195.2069.peg.1050... *
        >Cluster 490
        0	551aa, >fig|195.2045.peg.1513... *
        1	530aa, >fig|195.2053.peg.1024... at 84.34%
        2	536aa, >fig|195.2186.peg.1791... at 94.40%
        3	541aa, >fig|195.2242.peg.678... at 94.82%
        >Cluster 491
        0	501aa, >fig|195.2049.peg.1588... at 99.80%
        1	551aa, >fig|195.2166.peg.1287... *
        2	453aa, >fig|195.2206.peg.1626... at 99.56%
        3	453aa, >fig|195.2226.peg.1423... at 99.78%

        Each cluster has a number, here representing a CDS; each entry starts with a >Cluster line followed by the cluster number
        Under each cluster there are all entries for that cluster in each sample genome, based on a sequence identity (here 80% from pipeline)
        One entry under a cluster is the representative sequence, marked with an asterisk (*), this is the entry we want to get its gene name from PATRIC 
        p.s., 195.2029.peg.1780 is the gene name for cluster 0

        The cluster parser aims to match each cluster number to its representative sequence and get the gene name from PATRIC using api requests
    '''
    with open(clstr_path, 'r') as clstr_file:

        found_representative=False #flag to know if we found a representative for the clstr number
        cluster='' #initialize cluster number

        cluster_dict={} #initialize dictionary to store cluster number and representative sequence

        for line in clstr_file:

            if line.startswith(">"): #if its a cluster number file, save the cluster number for the next iteration
                cluster=line[1:].strip()
                found_representative=False #reset to false in order to find rep

            elif not found_representative:
                if "*" in line:
                    found_representative=True
                    representative=line.split(",")[1].replace('... *', '').strip()[1:]
                    representative=representative.split("|")[1]
                    # print(representative)
                    cluster_dict[cluster]=representative

            else:
                continue

    df=pd.DataFrame(cluster_dict.items(), columns=['cluster', 'gene_representative'])
    df.set_index('cluster', inplace=True)
    df.index.name='Cluster'
    
    return df

def get_representative_products(clstr_fasta_path:str)->pd.DataFrame:
    '''
    Takes a CD-HIT fasta output and returns a dataframe that matches each CDS PATRIC ID with its product name

    the fasta output is of this form:

    >fig|195.2024.peg.83 Putative oxidoreductase ferredoxin-type protein, clusters with CPO
    MNFSQISDACVKCGKCIPVCTIHEVNRDETTSPRGFLDLLAAYKEEKLELDKEAKKIFES
    CFLCTNCVEVCPSKLRVDNVIEEVRYDIAKKFGIAWYKKIIFFFLRRRKILDLVAKLGYV
    FQSCAFKIQSQDQNVGMKAKFSMPFVKKGRLLTSFNKKSFLNSNPDFIDNGGEKTVGFFV
    GCLANYFYIDTANAVLKIAKEVKINVDLMKEQVCCGAPQFFTGDFKSVEILAKKNIEYFE
    KKLEKLDAIIIPEATCSAMLKIDLEHFFNMQNEPEWAKRAQKISSRIYMASEYFYKFTNL
    KELLESKKKLNYSITYHDPCHARKMQGVFKEPRELLKANYHFVEMSNPNACCGFGGVSMQ
    TDYYDRALSVGLKKASMIDESKACVVSAECSACRMQISNALEQNSSKAIFASPLELIAKA
    L
    >fig|195.2024.peg.542 Septum-associated rare lipoprotein A
    MKPYTINGKTYYPTVVSVGETADGIASWYGPGFHGKKTSNGETYNQNGLTAAHKTLPMNT
    ILKVTNLNNNRQVTVRVNDRGPFVNNRIIDLSKGAASQIDMIAAGTAPVRLEVIGFGSAN
    SGNNVVHSNINYGASGGIANNGQIYEGGNFMVQIGAFKNPSGAQTIASRYKTYRTYSSTI
    RKSSVDGLSRVFLTGFRSEEEARDFAASGAFAGAFVVRE

    The aim is to match 195.2024.peg.542 with Septum-associated rare lipoprotein A in teh dataframe

    param:
    ------
    - clstr_fasta_path: str, path of the fasta output

    return:
    -------
    - df: pd.DataFrame, df of columns gene_representative and product_name
    '''

    dict={}

    with open(clstr_fasta_path,'r') as f:
        for line in f.readlines():
            if line.startswith(">"):
                pattern=">fig.(\d+\.\d+.peg.\d+) (.+)"
                gene_representative=re.match(string=line, pattern=pattern).group(1)
                product_name=re.match(string=line, pattern=pattern).group(2)
                dict[gene_representative]=product_name
    

    df = pd.DataFrame(dict.items(), columns=['gene_representative', 'product_name'])

    # -- autoincrement the product name, accounting for alleles
    autoincrement_col(df, 'product_name')
    
    return df

def combine_cluster_product(clstr_rep_df:pd.DataFrame, rep_prod_df:pd.DataFrame)->pd.DataFrame:
    '''
    Combines the cluster representative dataframe with the product name dataframe

    param:
    ------
    - clstr_rep_df: pd.DataFrame, dataframe with cluster number and gene representative (get_cluster_representatives() output)
    - rep_prod_df: pd.DataFrame, dataframe with gene representative and product name (get_representative_products() output)

    return:
    -------
    - df: pd.DataFrame, dataframe with cluster number, gene representative and product name
    '''
    # if clstr rep df doesnt have a column named Cluster, create one out of index
    if 'Cluster' not in clstr_rep_df.columns:
        clstr_rep_df['Cluster']=clstr_rep_df.index

    df=pd.merge(clstr_rep_df, rep_prod_df, on='gene_representative', how='left')
    df.drop(columns='gene_representative', inplace=True)

    df.set_index('Cluster', inplace=True)

    return df

def get_cluster_pan_gene_class(pan_annot_path:str)->pd.DataFrame:
    '''
    Takes a pangenome annotation file and returns a dataframe with cluster number and gene class
    This annotation file is a <species>_pangenome.csv file output from the pangenome analysis pipeline

    param:
    ------
    - pan_annot_path: str, path to the pangenome annotation file

    return:
    -------
    - pan_df: pd.DataFrame, dataframe with cluster number and gene class
    '''
    with open(pan_annot_path, 'r') as f:
        pan_df=pd.read_csv(f, index_col=0)
        pan_df.drop(pan_df.columns[[1]], axis=1, inplace=True)
        pan_df.columns=['Cluster', 'pan_gene_class']

    pan_df.set_index('Cluster', inplace=True)
    return pan_df

def get_cluster_representatives(clstr_path:str, save_csv:bool=False)->pd.DataFrame:
    '''
    Take a CD-HIT clstr output file and return a dictionary with cluster number as key and representative sequence as value.
    Optionally saves the dictionary as a csv file.

    parameters:
    ----------
    clstr_path: str, path to CD-HIT clstr file
    save_csv: bool, if True saves the dictionary as a csv file (optional)

    return:
    -------
    df: pd.DataFrame, data frame with cluster number as key and representative sequence as value

    description:
    ------------

    CD-HIT cluster file parser:
        .clstr file format is:

        >Cluster 489
        0	552aa, >fig|195.2069.peg.1050... *
        >Cluster 490
        0	551aa, >fig|195.2045.peg.1513... *
        1	530aa, >fig|195.2053.peg.1024... at 84.34%
        2	536aa, >fig|195.2186.peg.1791... at 94.40%
        3	541aa, >fig|195.2242.peg.678... at 94.82%
        >Cluster 491
        0	501aa, >fig|195.2049.peg.1588... at 99.80%
        1	551aa, >fig|195.2166.peg.1287... *
        2	453aa, >fig|195.2206.peg.1626... at 99.56%
        3	453aa, >fig|195.2226.peg.1423... at 99.78%

        Each cluster has a number, here representing a CDS; each entry starts with a >Cluster line followed by the cluster number
        Under each cluster there are all entries for that cluster in each sample genome, based on a sequence identity (here 80% from pipeline)
        One entry under a cluster is the representative sequence, marked with an asterisk (*), this is the entry we want to get its gene name from PATRIC 
        p.s., 195.2029.peg.1780 is the gene name for cluster 0

        The cluster parser aims to match each cluster number to its representative sequence and get the gene name from PATRIC using api requests
    '''
    with open(clstr_path, 'r') as clstr_file:

        found_representative=False #flag to know if we found a representative for the clstr number
        cluster='' #initialize cluster number

        cluster_dict={} #initialize dictionary to store cluster number and representative sequence

        for line in clstr_file:

            if line.startswith(">"): #if its a cluster number file, save the cluster number for the next iteration
                cluster=line[1:].strip()
                found_representative=False #reset to false in order to find rep

            elif not found_representative:
                if "*" in line:
                    found_representative=True
                    representative=line.split(",")[1].replace('... *', '').strip()[1:]
                    representative=representative.split("|")[1]
                    # print(representative)
                    cluster_dict[cluster]=representative

            else:
                continue

    df=pd.DataFrame(cluster_dict.items(), columns=['cluster', 'gene_representative'])
    df.set_index('cluster', inplace=True)
    df.index.name='Cluster'

    return df

# -------------- extra extra attributes (df) --------------------------

clstr_patric_id_df=get_cluster_representatives(clstr_path)
clstr_gene_class_df= get_cluster_pan_gene_class(pan_annot_path)
_patric_id_product_df=get_representative_products(clstr_fasta_path)
clstr_product_df=combine_cluster_product(clstr_patric_id_df, _patric_id_product_df)
clstr_freq_df = get_cluster_frequency(freq_path)

# -- getters 

def get_clstr_patric_id_df():
    return clstr_patric_id_df

def get_clstr_gene_class_df():
    return clstr_gene_class_df

def get_patric_id_product_df():
    return _patric_id_product_df

def get_clstr_product_df():
    return clstr_product_df

def get_clstr_freq_df():
    return clstr_freq_df

# -----------------------------------------------------------------------

def get_cluster_attributes(clstr_path=clstr_path, freq_path=freq_path, pan_annot_path=pan_annot_path, clstr_fasta_path=clstr_fasta_path):

    clstr_patric_id_df=get_cluster_representatives(clstr_path)
    clstr_gene_class_df= get_cluster_pan_gene_class(pan_annot_path)
    _patric_id_product_df=get_representative_products(clstr_fasta_path)
    clstr_product_df=combine_cluster_product(clstr_patric_id_df, _patric_id_product_df)
    clstr_freq_df=get_cluster_frequency(freq_path)

    df=pd.concat([clstr_patric_id_df, clstr_gene_class_df, clstr_product_df, clstr_freq_df], axis=1)

    #if there exist a Cluster col remove
    if 'Cluster' in df.columns:
        df.drop('Cluster', axis=1, inplace=True)
    
    return df

def create_dicts(df:pd.DataFrame):
    '''
    takes a nxm data frame and returns m dictionaries of keys the index, and values each col repectively

    param:
    -----
    - df: pd.DataFrame, index=Cluster, cols=gene_representative,	pan_gene_class,	product_name,	frequency

    returns:
    --------
    - dicts: dict, dict of m dictionaries - key is the col name, value is the dictionary of the col
    '''
    dicts = {col: df[col].to_dict() for col in df.columns}
    return dicts

def set_node_attributes_by_cluster(df:pd.DataFrame, G:nx.graph):
    '''
    takes a dataframe and a graph and sets the node attributes for each cluster

    param:
    -----
    - df: pd.DataFrame, index=Cluster, cols=gene_representative,	pan_gene_class,	product_name,	frequency
    - G: nx.Graph, graph to set the node attributes

    returns:
    --------
    - G: nx.Graph, graph with node attributes set
    '''
    dicts = create_dicts(df)

    for key in dicts.keys():
        nx.set_node_attributes(G, dicts[key], key)

    return G





def set_cluster_attributes(G:nx.Graph, pipeline_output=_output_path)->nx.Graph:
    '''
    takes a Graph and the pipeline's output path and performs all of the following: 

    - get all the files paths required for the attributes
    - generate df for each and concatenate them all in one of index Cluster
    - create dictionaries out of it and perform nx.set_node_attributes using each on G

    param:
    ------
    G: nx.Graph, network of nodes Clusters

    return:
    -------
    G: nx.Graph, network of nodes Clusters, attributes: pan_gene_class, gene_representative, pan_gene_class, frequency
    '''

    p1, p2, p3, p4 = generate_files_paths(pipeline_output)
    df = get_cluster_attributes(p1, p2, p3, p4)
    G = set_node_attributes_by_cluster(df, G)

    return G

# --------------- LOR calculation ----------------

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
    pheno_df = pd.read_csv(pheno_path, index_col=0, dtype={0: str})

    # -- already processed, this will give errors
    # if presence_df is not None:
    #     samples_presence= presence_df.columns
    #     mask = pheno_df.index.astype(str).isin(samples_presence)
    #     pheno_df = pheno_df.loc[mask]

    return pheno_df

def split_matrix_by_phenotype(unlabeled_presence_df:pd.DataFrame, pheno_df: pd.DataFrame)->(pd.DataFrame, pd.DataFrame):
    '''
    takes a gene absence presence dataframe (n x m) with the samples classification and splits it to 2 dataframes:
        - one that represents all R samples (samples found in R_list)
        - one that represents all S samples (found in S_list)
    and returns them in this order: R then S

    param:
    ------
        - unlabeled_presence_df: (pd.DataFrame) the data frame (GxS)
        - R_list: (list) samples genome_id that have a R (1) phenotype
        - S_list: (list) samples having the 0 phenotype
    
    return:
    ------
        - R_df: (pd.DataFrame) data frame where columns are only for R samples
        - S_df: (pd.DataFrame) df where cols are only for S samples

    '''

# -- already taken care of
    # # make sure pheno_df doesnthave more samples than those specified in the presence_df, subsample it otherwise
    # samples_presence= unlabeled_presence_df.columns
    # mask = pheno_df.index.astype(str).isin(samples_presence)
    # pheno_df = pheno_df.loc[mask]

    # unlabeled_presence_df.columns=unlabeled_presence_df.columns.astype('float')
    # pheno_df.index=pheno_df.index.astype('str')

    # print(pheno_df.index)
    # print(unlabeled_presence_df.columns)


    # get the list of R and S samples:
    R=[];S=[];U=[]
    for sample in pheno_df.index:

        if pheno_df.loc[sample].values[0]==1:
            # sample=str(sample)
            R.append(sample)
        elif pheno_df.loc[sample].values[0]==0:
            # sample=str(sample)
            S.append(sample)
        else:
            sample=str(sample)
            U.append(sample)

    # All we care for care R and S which designate the list of resistant and susceptible samples respectively

    # print(R)

    # print(R[0])
    # print('562.22603' in unlabeled_presence_df.columns)

    R_df = unlabeled_presence_df[R]
    S_df = unlabeled_presence_df[S]

    # print(unlabeled_presence_df)

    return R_df, S_df

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

    indices=R.index

    new_df=pd.DataFrame({'R_present':R_present, 'R_absent':R_absent, 'S_present':S_present, 'S_absent':S_absent})
    new_df.index=indices

    #the 0.5 correction:
    for index in new_df.index:
        row=list(new_df.loc[index])
        if 0 in row: #check if any of the cols have value 0
            new_df.loc[index]=new_df.loc[index]+0.5

    # print(new_df)

    return new_df

def compute_log_odds_ratio(RS_counts_df:pd.DataFrame)->pd.DataFrame:
    '''
    takes the RS_counts_df and computes the log odds ratio for each gene in the dataframe

    param:
    ------
        - RS_counts_df: pd.DataFrame


    return:
    ------
        - log_odds_df: pd.DataFrame
    '''

    # n_R=len(R); n_S=len(S)


    R_present=RS_counts_df['R_present']
    R_absent=RS_counts_df['R_absent']
    S_present=RS_counts_df['S_present']
    S_absent=RS_counts_df['S_absent']

    log_odds=np.log((R_present/R_absent)/(S_present/S_absent))

    # print(log_odds)

    df=pd.DataFrame({"log_odds":log_odds})
    return df

def get_cluster_resistance_LOR(presence_df:pd.DataFrame, pheno_path:str)-> pd.DataFrame:
    '''
    Knits all the steps from splitting the matrix to computing the log odds ratio for each gene in the matrix

    param:
    -----
        * presence_df: pd.DataFrame, GxS
        * pheno_df: str, path to the phenotype file

    return:
    -------
        * log_odds_df: pd.DataFrame, Gx1
    '''
    pheno_df = generate_phenotype_df(pheno_path, presence_df)
    # print(pheno_df)
    
    R,S=split_matrix_by_phenotype(presence_df, pheno_df)
    RS_counts = generate_RS_presence_counts(R, S)
    log_odds_df = compute_log_odds_ratio(RS_counts)

    return log_odds_df

def transform_cluster_to_product(cluster_name, clstr_product_df=clstr_product_df):
    '''
    takes a cluster name and returns the product name of the cluster

    param:
    ------
    - cluster_name: str, cluster name
    - clstr_product_df: pd.DataFrame, df with cluster number and product name

    return:
    -------
    - product: str, product name of the cluster
    '''
    #match it ffrom clstr_poduct_df from the index
    cluster = clstr_product_df.loc[cluster_name]
    product=cluster['product_name']
    
    return product

def transform_product_to_cluster(product_name, clstr_product_df=clstr_product_df):
    '''
    takes a product and returns the cluster

    param:
    ------
    - product_name: str, product name
    - clstr_product_df: pd.DataFrame, df with cluster number and product name

    return:
    -------
    - cluster: str, product name of the cluster
    '''
    #match it ffrom clstr_poduct_df from the index
    product = clstr_product_df.loc[product_name]
    cluster=cluster['Cluster']
    
    return product

