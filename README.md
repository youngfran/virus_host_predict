# Virus Host Predict 
This repository contains code and data  to predict host taxonomic information from viral genomes.  

## Data  
Data was down loaded from https://www.genome.jp/virushostdb/ on the 25_1_2019.  
Virushostdb.tsv > inputs/VHDB_25_1_2019.csv 

The data dierectory contains 3  sub direcories.   
1. Nucleic acid. Refseq for each virus    > data/fna  
2. Amino acid sequence for the CDS regions. For each virus > data/faa  
3. Pfam domains - extract from the above aa sequences with  HMMER > data/pfs  

## Inputs  

CSV files containing the information for each dataset run with the corresponding results are in the results folder.  

## Code  
Contains the notebooks used to train and test the different datasets.  
The  ‘run_BACTt_DNA.ipynb notebook contains the best annotations.  
Each notebook generates one csv file in the results folder.  
