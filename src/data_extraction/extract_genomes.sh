#!/bin/sh

# if [ ! -f "data/genome_summary.csv" ]; then
#     curl ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_summary -o data/genome_summary.tsv
#     sed 's/\t/,/g' data/genome_summary.csv > data/genome_summary.csv
#     rm -f data/genome_summary.tsv
# fi  

# python src/data_Extraction/extract_accessions.py 

for file in `ls data/PATRIC_IDs/`;do

    dir_name=$(basename -s .txt $file)
    dir_path=data/genomes/$dir_name/

    if [ ! -d $dir_path ]; then
        mkdir $dir_path
    fi

    log=$dir_path/output.log
    echo "- Starting genomes from $file accession IDs" >> $log

    for i in `cat data/PATRIC_IDs/$file`; do
        if [ ! -f ${dir_path}${i}.fna ]; then 
            curl -s -o ${dir_path}${i}.fna "ftp://ftp.bvbrc.org/genomes/$i/$i.fna";
            echo "  - Downloaded $i" >> $log
        else
            echo "  - $i already exists" >> $log
        fi
    done

    echo "- Finished genomes from $file accession IDs" >> $log
done

echo "" >> $log
echo "---------------------------------------------------------------------------------------------------------" >> $log
echo "" >> $log