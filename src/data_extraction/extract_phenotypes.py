'''
This script downloads the phenotypes xlsx file from this paper's supplementary info

Aim is to create a space in data/phentoypes for each species
& save the SIR readings of each species-drug combination in a csv file
The csv file will have as rownames the genome ids and only one column which is the AMR SIR reading 
SIR will be further processed in the next step into a binary reading

------------------------------------------
'''

import os
import sys
import pandas as pd
import requests

def binarize(cell):

    zero=["susceptible", "non_resistant", "Susceptible-dose dependent", "susceptible*", "non_resistant*", "Susceptible-dose dependent*"]
    one= ["intermediate","resistant","non_susceptible","IS","intermediate*","resistant*","non_susceptible*","IS*"]
    
    if str(cell) in zero:
        return 0
    elif str(cell) in one:
        return 1
    else:
        return cell


def get_AMR_df()->pd.DataFrame:
    '''
    Downloads the xlsx file from the supplementary info 2 and returns the dataframe
    the values will be binarized based on the legend in the paper's xlsx file:
    "susceptible" "non resistant" and "Susceptible-dose dependent" will be 0, others 1
    '''
    if not os.path.exists('data/phenotypes'):
        os.makedirs('data/phenotypes')

    #download the xlsx file from the supplementary info 2
    file_path='data/SIR_readings.xlsx'
    link="https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-43549-9/MediaObjects/41467_2023_43549_MOESM5_ESM.xlsx"
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        r = requests.get(link)
        with open(file_path,'wb') as f:
            f.write(r.content)
    df = pd.read_excel(file_path, sheet_name='Consolidated_SIRs')
    #set the genome_id col as rownames:
    df.set_index('genome_id', inplace=True)

    return df

def create_species_drugs_files(df:pd.DataFrame):
    '''
    Takes a dataframe of cols: genome_id, species, and several drugs
    filter them to put them in a file for each species, each drug named species_drug.csv
    of this form: genome_id, SIR
    
    --------
    param:
        - df: the dataframe of the AMR SIR readings
    return:
        None
    '''
    species=df['species'].unique()
    drugs=list(df.columns[1:])
    for s in species:

        for d in drugs:
            sp_df=df[df['species']==s] #filtering a new df for only this species
            sp_df=sp_df[[d]] #drop all columns except the col d
            sp_df=sp_df.dropna(subset=[d]) #drop all NaN in col d
            
            if sp_df.empty:
                continue
            
            sp_df.columns=['SIR']
            sp_df=sp_df['SIR'].apply(binarize)

            sp_df.to_csv(f'data/phenotypes/{s.replace(" ","_")}_{d}.csv')
    return None


def main():
    df=get_AMR_df()
    create_species_drugs_files(df)

if __name__=='__main__':
    main()

