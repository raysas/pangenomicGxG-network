#!/bin/bash

# * streptomycin
# * sulfamethoxazole
# * tetracycline
# * cefalothin
# * trimethoprim
# * amoxicillin
# * ampicillin
# * levofloxacin
# * ciprofloxacin

species='Escherichia_coli'

for antibiotic in streptomycin sulfamethoxazole cefalothin tetracycline trimethoprim amoxicillin ampicillin levofloxacin ciprofloxacin; do
    mkdir ../../results/${species}_${antibiotic}
done

# -- cpying the amoxicillin analysis notebook and replacing the antibiotic name
for antibiotic in streptomycin sulfamethoxazole cefalothin trimethoprim tetracycline ampicillin levofloxacin ciprofloxacin; do
    cp Escherichia_coli_amoxicillin_analysis.ipynb Escherichia_coli_${antibiotic}_analysis.ipynb
    sed -i "s/amoxicillin/${antibiotic}/g" Escherichia_coli_${antibiotic}_analysis.ipynb
done

# -- rm if needed to cp again
for antibiotic in streptomycin sulfamethoxazole cefalothin trimethoprim amoxicillin ampicillin levofloxacin ciprofloxacin; do
    rm Escherichia_coli_${antibiotic}_analysis.ipynb
done