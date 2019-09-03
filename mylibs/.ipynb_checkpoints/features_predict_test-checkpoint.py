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

from collections import Counter
import seaborn as sns


def get_class_lists(subset,viruses,hosts):
# returns two list of viruses one for each  positive class and negative class    
    (label,label_tax,pool,pool_tax,baltimore) = subset
    if baltimore =='all':
        baltimore = 'NA'  # all classes (except satallites)
    
    # Get a list of all the viruses in the labelled class and the rest of the pool
    bc = [[],[]]
    for v, vd in viruses.items() :
        if baltimore in vd['class']:
            host_labels = [hosts[h][label_tax] for h in vd['hosts'] if hosts[h][pool_tax] == pool]    
            if label in host_labels:
                bc[0].append ((v,label))
            else:
                bc[1].append((v,'Other'))                     

    # Get a random sample for each class of size n , the size of smallest class                    
    n = min(len(bc[0]),len(bc[1]))
    n=10 #  for testing
    datas =  []
    for clss in bc:
        data = (random.sample(clss, n))
        datas.append(data)
    
    return datas,n







def split_data(viruses,data,n):

    if n < 50:
        split = 0.75
    else:
        split = 0.8 # ratio of training:test (traininig includes validatio data)

    datasets = {'training':{},'test':{}}
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

    for v,l in trainingSet:
        datasets ['training'][v]= {'label':l, 'refseqs':viruses[v]['refseqs']}
    for v,l in testSet:
        datasets ['test'][v]= {'label':l , 'refseqs':viruses[v]['refseqs']}
    print (f' size of training and test set {len(trainingSet)},  {len(testSet)}')
    return (datasets)
    


# ## Functions to create feature matrix


def SymbolToNumber (symbol,symbols):
    return symbols[symbol]


def PatternToNumber (pattern,k,symbol_dict):
    mod = symbol_dict['mod']
    number = 0
    for i in range (0,k):
        n = SymbolToNumber ( pattern [k-1-i],symbol_dict)
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

# 
def compliment(dna):
    """Return the reverse compliment of a sequence """
    bp = { 'a':'t', 't':'a', 'c':'g','g':'c'}
    revdna =''    
    for i in range(len(dna)-1,-1,-1):
         revdna += bp[dna[i]]
    return (revdna)
#

    
    
    
def get_word_frequ(v,refs,k,vocab_len,symbol_dict, faapath):
    """returns the word frequency array for all genomes for a virus """
    word_freqs = np.zeros(shape =(vocab_len),dtype=np.int32)
    
    for ref in refs:
        
        fastaname = os.path.join(faapath, ref.strip()+".fasta")
 
        try:
            with open(fastaname, 'rt') as handle:

                for record in SeqIO.parse(handle, "fasta"):
                    sequence = record.seq 
                    for i in range(len(sequence)-k+1):
                        word = sequence[i:i+k]
                        if check_word(word,symbol_dict):
                            kmer_index = PatternToNumber (word,k, symbol_dict)
                            #print (word, k,kmer_index)
                            word_freqs [ kmer_index] += 1
                            # If this for DNA/RNA then add the k-mers complinent
                            if symbol_dict['mod'] == 4:
                                comp_kmer_index = PatternToNumber(compliment(word),k,symbol_dict)
                                word_freqs [ comp_kmer_index] += 1
        except IOError:
            #print('input error')
            pass
        
    return word_freqs



def DNA_features (datasets,k,symbol_dict,filepath):
    Xs,v_index,f_index = AA_features (datasets,k,symbol_dict,filepath)
    return tuple(Xs),v_index,f_index


def AA_features (datasets,k,symbol_dict,filepath):
    vocab_len =  symbol_dict['mod'] **k
    empty = []
    f_index = get_feature_names(symbol_dict,k)
    
    Xs = [None]*2
    v_index =[None] *2
    for index_ds,(key,dataset) in enumerate (datasets.items()):
        
        v_index [index_ds] = []
        virus_index = 0
        X = np.zeros((len(dataset) ,vocab_len))
        print (f'start of making  {key} X,{len(dataset)}, {np.shape(X)}')
        for v,vdict in dataset.items():
            word_freqs = get_word_frequ( v, vdict ['refseqs'],k, vocab_len,symbol_dict, filepath)
            if (np.sum(word_freqs)) != 0:
                for i in range (vocab_len): 
                    X [virus_index,i] = word_freqs[i]

                v_index [index_ds].append(v)
                virus_index += 1
            
        
        #remove empty rows from missing viruses
        X1 = X[:len(v_index[index_ds]),:]
        # normalise for length of genome(s)
#        XN = np.divide(X1,X1.sum(axis = 1)[:,None])
        print (len (v_index[index_ds]),np.shape(XN))
                    
        Xs[index_ds] = X1
        
    
    return Xs,v_index,f_index
    


# In[9]:


def PC_features (datasets,k,symbol_dict,filepath):
    Xs,v_index = AA_features (datasets,k,symbol_dict,filepath)
    f_index = get_feature_names(symbol_dict,k)
    return tuple(Xs),v_index,f_index


# In[10]:


def PFAM_list(datasets,dom_dict):
# get unique list o fPFAM ids for feature matrix

    PFAMs = []
    for key,dataset in datasets.items():
        for v in dataset:
            doms = dom_dict[v]
            for d in doms:
                PFAMs.append(d) 

    PFams = list(set (PFAMs)) 
    return PFams


# In[11]:


def make_PFAM_dict(datasets,filepath):
    """ Create a PFAM dictionary for all the viruses in both training and test datasets. 
         Returns {taxid:[PFAM1,..],..}
         Input, dataset dictionaries ,filepath where pfs file are stored
         output domain_dict"""


    domain_dict={}
    count = 0
    missing =0
    
    for key,dataset in datasets.items():
        for  virus,v_dict in dataset.items():
            for ref in v_dict['refseqs']:
                filename = os.path.join(filepath, ref.strip()+".pfs")
                
                if not os.path.isfile(filename): 
                    #print (filename,'Isfile:file not found')
                    missing += 1
                else:  
                    doms =[]
                    with open(filename, "r") as f:
                        line = f.read()
                        for d in line.split():
                            doms.append(d)

                    if virus in domain_dict:
                        domain_dict[virus].extend(doms)
                        # print( virus,'multiple refs')
                    else:
                        domain_dict[virus] = doms
                        count += 1
                 
    print ('There are {} viruses with domains, {} viruses missing domain info'.format(count,missing),  )
    return domain_dict


def get_domain_features(datasets,dom_dict):
# remove viruses with no domain info from dataset used

    ds= {}
    
    for key,dataset in datasets.items():
        ds[key] ={}
        for virus,vdict in dataset.items():
            if virus in dom_dict:
                ds[key][virus]=vdict
    
         
    PFams = PFAM_list(ds,dom_dict)
    
    Xs = [None]*2
    v_index =[None] *2
    for index_ds,(key,dataset) in enumerate (ds.items()):
        v_index [index_ds] = [None] * len(dataset)
        X = np.zeros((len(dataset) ,len(PFams)))
        for i,(v,vdict) in enumerate(dataset.items()):
            v_index[index_ds][i] = v
            doms = dom_dict[v]
            for dom in doms:
                j = PFams.index(dom)
                X[i,j]+= 1 
        XN = np.divide(X,X.sum(axis = 1)[:,None]) 
        print (len(dataset), np.shape(XN))
        Xs[index_ds] = XN
        
    
    return Xs,v_index,PFams
    

def Domain_features(dataset,bla_1,bla_2,filepath):
    domain_dict = make_PFAM_dict(dataset,filepath)
    Xs,v_index,PFams = get_domain_features(dataset,domain_dict)
    return tuple(Xs) ,v_index, PFams


# ### Get the labels y
# In binary (True or false). For the virus index list (missing viruses removed

# In[14]:


def  get_labels(label,v_index,datasets):
    Ys = []
    for index_ds,(key,dataset) in enumerate (datasets.items()):
        
        y = [None] * len(v_index[index_ds]) 
        for i,v in enumerate(v_index[index_ds]):
            if  dataset[v]['label'] ==  label:
                #print ('true')
                y[i] = 1
            else:
                y[i] = 0
        Ys.append(y)
        print ( f'labels y {index_ds},{len(y)}')
    
    return tuple(Ys)


def test_prediction(X_train,y_train,X_test,y_test):
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
        with open(csvfile, 'a') as csvfile:
            fieldnames = ['positive label','label tax group','pool label','pool tax group','Baltimore', 'N in class' , 'Features','k','AUC',]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(results)  
        
  