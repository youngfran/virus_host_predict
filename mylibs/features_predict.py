 # coding: utf-8


import pickle
import numpy as np
import pandas as pd
import os
import csv
import string
import random 

from Bio import SeqIO


from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score,roc_auc_score
from sklearn.metrics import confusion_matrix

from collections import Counter, defaultdict

   
def get_class_lists(subset,viruses,hosts):
# returns two list of viruses one for each  positive class and negative class    
    (label,label_tax,pool,pool_tax,baltimore) = subset
    if baltimore =='all':
        baltimore = 'NA'  # all classes (except satallites)
    
    # Get a list of all the viruses in the labelled class and the rest of the pool
    both_classes = [[],[]]
    for v, vd in viruses.items() :
        if baltimore in vd['class']:
            host_labels = [hosts[h][label_tax] for h in vd['hosts'] if hosts[h][pool_tax] == pool] 
            if len(host_labels) > 0:
                if label in host_labels:
                    both_classes[0].append ((v,label))
                else:
                    both_classes[1].append((v,'Other'))                     

    # Get a random sample for each class of size n , the size of smallest class                    
    
    n = min(len(both_classes[0]),len(both_classes[1]),nmax)
#    n=25 #  for testing
    datas =  []
    for clss in both_classes:
        data = (random.sample(clss, n))
        datas.append(data)
    
    return datas,n



def split_data(viruses,data,n):

    if n < 50:
        split = 0.75
    else:
        split = 0.8 # ratio of training:test (traininig includes validatio data)

    
    trainingSet =  []
    testSet =[]
    random.seed(10)
    
    for viruslist in data:
        for i in range(n):
            if random.random() < split:
                trainingSet.append(viruslist[i]) 
            else:
                testSet.append(viruslist[i]) 
    random.shuffle(trainingSet)
    random.shuffle(testSet)
    
    datasets = {'training':{},'test':{}}
    for v,l in trainingSet:
        datasets ['training'][v]= {'label':l, 'refseqs':viruses[v]['refseqs']}
    for v,l in testSet:
        datasets ['test'][v]= {'label':l , 'refseqs':viruses[v]['refseqs']}
    print (f' size of training {len(trainingSet)} and test set   {len(testSet)}')
    return (datasets)
    


# ## Functions to create feature matrix

def extract_kmers (sequences,k,symbol_dict):
    if k == 0: # domainsdom_list = list({dom for s in sequence.values() for dom in s})
        
        X,seq_index,f_index = get_domains(sequences,symbol_dict)
    else:  # everything else
        vocab_len =  symbol_dict['mod'] **k
        f_index = get_feature_names(symbol_dict,k)
        seq_index = []
        X = np.zeros((len(sequences) ,vocab_len))

        for i,(accn, seq) in enumerate (sequences.items()):
            X[i] = get_word_frequ (seq,k, vocab_len,symbol_dict)
            seq_index.append(accn)
# normalise for length of genome(s)
    XN = np.divide(X,X.sum(axis = 1)[:,None])
    #print (np.shape(XN)) 
    return XN,seq_index,f_index

def get_domains(sequence,dom_list):
    ''' input: sequences dictionary of training and test sequence dicts'''
    seq_index = []
    X = np.zeros((len(sequence) ,len(dom_list)))
    for i,(accn, seq) in enumerate (sequence.items()):
        for dom in seq:
            j = dom_list.index(dom)
            X[i,j] += 1
        seq_index.append(accn)
    return X,seq_index,dom_list

def get_word_frequ(sequence,k,vocab_len,symbol_dict):
    """returns the word frequency array for all genomes for a virus """
    word_freqs = np.zeros(shape =(vocab_len),dtype=np.int32)
    for i in range(len(sequence)-k+1):
        word = sequence[i:i+k]
        if check_word(word,symbol_dict):
            kmer_index = patternToNumber (word,k, symbol_dict)
            word_freqs [ kmer_index] += 1
            # If this for DNA/RNA then add the k-mers complinent
            if symbol_dict['mod'] == 4:
                word_freqs [patternToNumber(compliment(word),k,symbol_dict)] +=1   
    return word_freqs

def get_sequences(faapath,datasets,ext):
    missing_files = []
    v_with_missing_files = []
    ext = f'.{ext}'
    sequences = {'training':{},'test':{}}
    for key,dataset in datasets.items():
        for v,vdict in dataset.items():
            sequence = defaultdict(str)
            for ref in vdict['refseqs']:
                filename = os.path.join(faapath, ref.strip()+ ext)
                try:
                    with open(filename, 'rt') as handle:
                        if ext == '.fasta':
                            for record in SeqIO.parse(handle, "fasta"):                   
                                sequence[v] += str(record.seq)
                        else: # pfs file - list of domains  
                            line = handle.read()
                            sequence[v] = [d for d in line.split()]
                                
                except IOError:
                    missing_files.append((v,ref))
                    
            if len(sequence) > 0:
                sequences[key].update(sequence)
            else:
                v_with_missing_files.append(v)
    print ( f'viruses number of missing files {len(set(v_with_missing_files))}')
    return (sequences)



def symbolToNumber (symbol,symbols):
    return symbols[symbol]

def patternToNumber (pattern,k,symbol_dict):
    mod = symbol_dict['mod']
    number = 0
    for i in range (0,k):
        n = symbolToNumber ( pattern [k-1-i],symbol_dict)
        number += n* mod**i 
    return number

def get_feature_names(symbols,k):
    all_kmers = []
    n2sym = { v:k for k,v in symbols.items()}
    mod = symbols['mod']
    for number in range(mod**k):
        pattern = ''
        for i in range (0,k):
            n = number% mod
            number = (number - n)/mod 
            symbol = n2sym[n]
            pattern = symbol + pattern
        all_kmers.append (pattern)
    return all_kmers

def check_word(mystring,symbols):
    """Check for illegal characters """
    return all(c in symbols for c in mystring)

def compliment(dna):
    """Return the reverse compliment of a sequence """
    bp = { 'a':'t', 't':'a', 'c':'g','g':'c'}
    revdna =''    
    for i in range(len(dna)-1,-1,-1):
         revdna += bp[dna[i]]
    return (revdna)

def  get_labels(label,v_index,dataset):  
    #print('y', label,len(v_index))
    for i,v in enumerate(v_index):
        y = [1 if dataset[v]['label'] ==  label else 0 for v in v_index ] 
        #print('y', label,len(v_index),y)
    return np.asarray (y)


def get_feature_matrices(sequences,datasets,label,k,symbol_dict):
   
    if k == 0: # domainsdom_list = list({dom for s in sequence.values() for dom in s})
        symbol_dict = list({dom for sequence in sequences.values() for s in sequence.values() for dom in s}) 
    X_train, seq_index, f_index = extract_kmers (sequences['training'],k,symbol_dict)
    Y_train =  get_labels (label,seq_index,datasets['training'])
 
    X_test, seq_index, f_index = extract_kmers (sequences['test'],k,symbol_dict)
    Y_test =  get_labels (label,seq_index,datasets['test'])
           
    return X_train,X_test,Y_train,Y_test

def test_prediction(X_train,X_test,y_train,y_test):
    clf = make_pipeline(StandardScaler(),svm.SVC(kernel='linear',probability=True) )
    clf.fit(X_train, y_train)
    y_pred_probs= clf.predict_proba(X_test)[:,1]
    AUC =round( roc_auc_score(y_test, y_pred_probs),3)        
    return {'AUC':AUC}


def results2CSV(results, subset, csvfile):
    (label,label_tax,pool,pool_tax,balt) = subset
    #results = { k : round(v, 3) for k,v in results.items()}
    results['positive label'] = label
    results['label tax group']= label_tax
    results['pool label']= pool
    results['pool tax group']= pool_tax
    results['Baltimore'] = balt
    
    if os.path.isfile(csvfile):
        with open(csvfile, 'a') as csvfile:
            fieldnames = ['positive label','label tax group','pool label','pool tax group','Baltimore', 'N in class' , 'Features','k','AUC' ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writerow(results) 
    else:
        with open(csvfile, 'w') as csvfile:
            fieldnames = ['positive label','label tax group','pool label','pool tax group','Baltimore', 'N in class' , 'Features','k','AUC',]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(results)  

