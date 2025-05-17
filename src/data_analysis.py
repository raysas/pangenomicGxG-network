'''
This one is to manipulate the different data frames that will be of need in several files.
Mainly:  
- getting presence matrices
- getting labeled matrices

Attributes 
-----------
- pan_path: the pangenome pipeline output path

Derived attributes:
- species: the species name
- parent_pan_path: the parent path of the pan_path
- genome_ids_file: the path to the genome ids file
- clstr_freq_file: the path to the cluster frequencies file
- clstr_file: the path to the cluster file
- clstr_fasta_file: the path to the cluster fasta file
- samples_df: the dataframe of the genome ids
- samples_list: the list of genome ids

'''

import os
import pandas as pd
import re 
import sys

os.chdir(os.path.expanduser('~/capstone-project'))
sys.path.append("src")

from cluster_analysis import *

os.chdir(os.path.expanduser('~/capstone-project'))

# ----------------------- Attributes -----------------------

pan_path = 'data/pangenome_pipeline_output/Escherichia_coli'
species = os.path.basename(pan_path)

parent_pan_path = os.path.dirname(pan_path) # used for genome_ids

# ----------------------- Getters and Setters -----------------------
def get_pan_path():
    return pan_path

def set_pan_path(path):
    global pan_path
    pan_path = path
    if pan_path[-1]=='/':
        pan_path=pan_path[:-1]

    if not os.path.exists(pan_path):
        raise FileNotFoundError(f"Pipeline output Path {pan_path} does not exist")
    
    global species
    species = os.basename(pan_path)

# ----------------------- Extra attributes -----------------------

def generate_presence_matrix(pan_path=pan_path, memoize=True):
    '''
    Generates a presence matrix from the cluster file and the genome ids file
    The presence file consensus file in this projet will be of the GxS format
    (Genes x Samples)  
    - rows: clusters (genes/alleles)
    - columns: genomes (samples)

    Parameters
    ----------
    pan_path: str
        The path to the pangenome pipeline output
    memoize: bool (optional)
        Whether to save the presence matrix to a csv file

    Returns
    -------
    pd.DataFrame
        The presence matrix
    '''

    species = os.basename(pan_path)
    parent_pan_path = os.path.dirname(pan_path)

    genome_ids_file=f"{parent_pan_path}/genome_ids/{species}_genome_ids.csv"
    clstr_freq_file=f"{pan_path}/{species}_cluster_frequencies.csv"
    clstr_file=f"{pan_path}/{species}.fasta.clstr"
    clstr_fasta_file=f"{pan_path}/{species}.fasta"

    samples_df = pd.read_csv(genome_ids_file, dtype=str)
    samples_list=samples_df['genome.genome_ids'].tolist()

    clstr_freq_df=pd.read_csv(clstr_freq_file, index_col=1)
    clstr_freq_df=clstr_freq_df.drop(clstr_freq_df.columns[0], axis=1)
    clstr_list=clstr_freq_df.index.tolist()

    df = pd.DataFrame(index=samples_list, columns=clstr_list)
    df = df.fillna(0)

    with open(clstr_file) as f:
        for line in f:
            if line.startswith(">Cluster"):
                #if its a > line get the cluster id and start counting its occurences in which samples by saving it to the local var
                cluster_id=line[1:].strip()
            else:
                #else its a line designating the sample genome, one that the last cluster_id is present in
                # get the 1st group matching here \d\t\d+aa, >fig\|([\d\.\d]+).+

                genome_id=re.match(r"\d+\t\d+aa, >fig\|([\d\.\d]+).+", line).group(1)
                genome_id=genome_id[:-1]  #removing the trailing .
                df.loc[genome_id, cluster_id]+=1

    df = df.T
    if memoize:
        df.to_csv(f'data/presence_matrices/{species}_GxS.csv')

    return df

def filter_presence_matrix(presence,pan_path=pan_path, unknowns=True, hypotheticals=True, freq_threshold=5, memoize=True):
    '''
    takes a presence matrix (df or path), a pangenome pipeline output path and filters out the clusters that are either:  
    - unknown
    - hypothetical
    - present in less than freq_threshold genomes
    All of these are optional

    returns a filetred SxG matrix (transposed - genes ready as features)

    Parameters
    ----------
    - presence: pd.DataFrame or str (GxS)
        The presence matrix or the path to the presence matrix
    - pan_path: str
        The path to the pangenome pipeline output
    - unknowns: bool (optional)
        True if we wanna filter them out
    - hypotheticals: bool (optional)
        True if we wanna filter them out
    - freq_threshold: int (optional)
        The frequency threshold to filter out clusters - 5 by default
    - memoize: bool (optional)
        Whether to save the filtered presence matrix to a csv file

    Returns
    -------
    - X_df: pd.DataFrame 
        The filtered presence matrix (SxG)
    - to_remove: list
        The list of clusters that were removed (maybe used to filter other things)
    '''

    clstr_df= get_cluster_representatives(f'{pan_path}/{species}.fasta.clstr')
    _product_df=get_representative_products(f'{pan_path}/{species}.fasta')
    product_df = combine_cluster_product(clstr_df, _product_df)
    freq_df = get_cluster_frequency(f'{pan_path}/{species}_cluster_frequencies.csv')

    # #  --counting those thar have hypotehtical in the product name
    product_df['hypothetical'] = product_df['product_name'].str.contains('hypothetical', case=False)
    # product_df['hypothetical'].sum()

    # #  --same for unknown
    product_df['unknown'] = product_df['product_name'].str.contains('unknown', case=False)
    # product_df['unknown'].sum()

    # now will get a list of clusters that are either hypo or unknown to filter out

    hypothetical_clusters=[] ; unknown_clusters=[]

    if hypotheticals:
        hypothetical_clusters = product_df[product_df['hypothetical']==True].index.tolist()
    if unknowns:
        unknown_clusters = product_df[product_df['unknown']==True].index.tolist()

    to_remove = list(set(hypothetical_clusters).union(unknown_clusters))

    # -- get the list of clusters that are present in less than 5 genomes
    low_freq_clusters = freq_df[freq_df['frequency']<freq_threshold].index.tolist()

    to_remove = list(set(to_remove).union(low_freq_clusters))


    if type(presence)==str:
        presence = pd.read_csv(presence, index_col=0) #if it's a path make it a df
    X_df=presence.T
    X_df.index = X_df.index.astype('float')

    X_df = X_df.drop(to_remove, axis=1)

    if memoize:
        X_df.to_csv(f'data/presence_matrices/{species}_filtered_SxG.csv')

    return X_df, to_remove

# def get_labeled_matrix(species=species, drug):
    