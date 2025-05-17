
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

presence_path="../../data/temp/X_Escherichia_coli_SxG_filtered_hypo_unknowns_freq5.csv" #add the method that generates presence to this script
presence_df= pd.read_csv(presence_path, index_col=0)
presence_df = presence_df.T #exceptionally because its a SxG
    # presence_df.index = presence_df.index.astype(str)



def generate_phenotype_df(pheno_path:str, presence_df:pd.DataFrame=None)->pd.DataFrame:
    '''
    takes in a path to a phenotype file to make a df out of it
    if a presence_df is provided, it will filter the pheno_df to only include samples that are in the presence_df

    param:
    ----
        - pheno_path: str
        - presence_df: pd.DataFrame (optional)

    return:
    ----
        - pheno_df: pd.DataFrame
    '''
    pheno_df = pd.read_csv(pheno_path, index_col=0)

    #samples stripped
    if presence_df is not None:
        samples_presence= presence_df.columns
        mask = pheno_df.index.astype(str).isin(samples_presence)
        pheno_df = pheno_df.loc[mask]

    return pheno_df

pheno_path="../../metadata/Escherichia_coli/Escherichia_coli_trimethoprim.csv"
pheno_df= generate_phenotype_df(pheno_path, presence_df)
    # pheno_df= pd.read_csv(pheno_path, index_col=0)


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
    print(indices)

    log_odds=np.log((R_copresent/R_coabsent)/(S_copresent/S_coabsent))

    print(log_odds)

    df=pd.DataFrame({"log_odds":log_odds}, index=indices)

    return df



def compute_cooccurence_LOR(gene_pairs:list, pheno_df:pd.DataFrame=pheno_df, presence_df:pd.DataFrame=presence_df )->pd.DataFrame:
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
    print('test1')
    R,S=split_samples_list_by_phenotype(pheno_df)
    print('test2')
    R_df,S_df=split_matrix_by_phenotype(presence_df, pheno_df)
    print('test3')
    R_cooccurence = get_cooccurence_matrix(R_df)
    print(R_cooccurence.shape)
    S_cooccurence = get_cooccurence_matrix(S_df)
    print('test5')
    numR=len(R); numS=len(S)

    print('split done and cooccurnece generated')
    print('-----------------------------------\n')

    RS_copresence_df=generate_RS_copresence_count(R_cooccurence, S_cooccurence, numR, numS, gene_pairs)
    print('copresence count done')
    print('-----------------------------------\n')

    log_odds_gene_pairs=compute_association_log_odds(RS_copresence_df)
    print('log odds done')
    print('-----------------------------------\n')

    return log_odds_gene_pairs

def main():
    # --
    # (f'../../data/temp/X_{species}_SxG_filtered_hypo_unknowns_freq5.csv', index_col=0)

    presence_path="../../data/temp/X_Escherichia_coli_SxG_filtered_hypo_unknowns_freq5.csv" #add the method that generates presence to this script
    presence_df= pd.read_csv(presence_path, index_col=0)
    presence_df = presence_df.T #exceptionally because its a SxG
    # presence_df.index = presence_df.index.astype(str)

    pheno_path="../../metadata/Escherichia_coli/Escherichia_coli_trimethoprim.csv"
    pheno_df= generate_phenotype_df(pheno_path, presence_df)
    # pheno_df= pd.read_csv(pheno_path, index_col=0)

    species = "Escherichia_coli"; drug="trimethoprim"

    gene_pairs=[]
    with open(f'../../data/temp/SVM_{species}_{drug}_top_1000_gene_pairs.txt', 'r') as f:
        for line in f:
            gene_pairs.append(tuple(line.strip().split("-")))

    gene_pairs=[("Cluster 0", "Cluster 20")]

    print('prepared gene pairs')
    print('-----------------------------------\n')

    LOR = compute_cooccurence_LOR(gene_pairs, pheno_df=pheno_df, presence_df=presence_df)

    print('saving LOR')
    print('-----------------------------------\n')

    LOR.to_csv(f'../../data/temp/LOR_SVM_{species}_{drug}.csv')

    print('fine')

if __name__ == "__main__":
    main()