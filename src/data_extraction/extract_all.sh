#!/bin/bash

#saving all the extracted data in subdirectories of data/
if [ ! -d data/ ]; then
    mkdir data/
fi

# 1. extract the genomes from PATRIC
## 1.1 download the Dataset1.json file from supplementary info
## 1.2 extract the genome IDs from the Dataset1.json file
## 1.3 download the genomes from the genome IDs
bash src/data_Extraction/PATRIC_extract_data.sh

# 2. create pangenomes for each species
## 2.1 run prokka to annotate the genomes - get GFF3 files
## 2.2 run roary to create the pangenome
## loop through the species
bash src/create_pangenome.sh data/genomes/Acinetobacter_baumannii/ Acinetobacter data/pangenomes/Acinetobacter_baumannii #running
bash src/create_pangenome.sh data/genomes/Campylobacter_coli/ Campylobacter data/pangenomes/Campylobacter_coli
bash src/create_pangenome.sh data/genomes/Campylobacter_jejuni/ Campylobacter data/pangenomes/Campylobacter_jejuni
bash src/create_pangenome.sh data/genomes/Enterobacter_cloacae/ Enterobacter data/pangenomes/Enterobacter_cloacae
bash src/create_pangenome.sh data/genomes/Enterococcus_faecium/ Enterococcus data/pangenomes/Enterococcus_faecium
bash src/create_pangenome.sh data/genomes/Escherichia_coli/ Escherichia data/pangenomes/Escherichia_coli
bash src/create_pangenome.sh data/genomes/Klebsiella_pneumoniae/ Klebsiella data/pangenomes/Klebsiella_pneumoniae
bash src/create_pangenome.sh data/genomes/Neisseria_gonorrhoeae/ Neisseria data/pangenomes/Neisseria_gonorrhoeae
bash src/create_pangenome.sh data/genomes/Pseudomonas_aeruginosa/ Pseudomonas data/pangenomes/Pseudomonas_aeruginosa
bash src/create_pangenome.sh data/genomes/Salmonella_enterica/ Salmonella data/pangenomes/Salmonella_enterica
bash src/create_pangenome.sh data/genomes/Staphylococcus_aureus/ Staphylococcus data/pangenomes/Staphylococcus_aureus
bash src/create_pangenome.sh data/genomes/Streptococcus_pneumoniae/ Streptococcus data/pangenomes/Streptococcus_pneumoniae

# 3. extract the phenotypes
python src/data_Extraction/extract_phenotypes.py