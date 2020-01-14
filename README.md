# Virus Host Predict 
This repository contains code and instructions to get the data  to predict host taxonomic information from viral genomes.  

## Data  
Data was downloaded from https://www.genome.jp/virushostdb/ on the 25_1_2019.    
Download the following two files:    
1. virushostdb.genomic.fna.gz  -> data/fna  
2. virushostdb.cds.faa.gz   -> data/fna 

Run the unpackfastafiles notebook to unpack the files. 
Unzip the pfs.zp in the data/pfs folder


The data directory contains 3  sub direcories.   
1. Nucleic acid. Refseq(s) for each virus    > data/fna  
2. Amino acid sequence for the CDS regions for each refseq > data/faa  
3. Pfam domains - extract from the above aa sequences with  HMMER  for each refseq> data/pfs  

## Inputs  

CSV files containing the information for each dataset run with the corresponding results are in the results folder.  
Virushostdb.tsv > inputs/VHDB_25_1_2019.csv 

## Code  
Contains the notebooks used to train and test the different datasets.  
The  run_TrainTest.ipynb notebook contains the best annotations. (select which datasets to run)  
Each notebook generates one csv file in the results folder.  

### Mylibs
1. vhdb.py - load virus host data base table as an object.
2. datasets.py  - class containing  dataset information.
3. featureset.py - class containing information for each featureset for the dataset object. 
