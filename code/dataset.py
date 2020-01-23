import os
import pandas as pd
import pickle
import random
import csv
from Bio import SeqIO
import sys
import numpy as np
from sklearn.model_selection import train_test_split 
from pathlib import Path
from featureset import FeatureSet
class DataSet (object):
    """ A dataset of positive and negative labelled viruses, with following attributes:
    label: string for name of positive class- host taxon group.
    tax: string with tax rank of label - one of ['phylum',.....,'species']
    pool: string with name of negative class
    pool_tax: tax of -ve class
    baltimore: string -baltimore class of viruses - one of ['all','RNA',DNA,(-)ssRNA,(+)ssRNA,dsRNA,
                                                            dsDNA,ssDNA,retro]
    featuresets: dict {fs_name: FeatureSet} fs_name eg. AA_4
    """ 
    
    # Class Attributes
    # Everthing in the same order as the features list
    # 
    feature_list = ['DNA','AA','PC','Domains']
    #filepaths = [Path('../data/fna/'),Path('../data/faa'),Path('../data/faa/'), Path('../data/pfs/')]
    #file_exts = ['.fasta','.fasta','.fasta','.pfs']
    filepaths = [(Path('../data/fna/'),'.fasta'),(Path('../data/faa/'),'.fasta'), \
                 (Path('../data/faa/'),'.fasta'), (Path('../data/pfs/'),'.pfs')]
    
    def __init__(self,subset,vh, feature_sets=None,tt_split=0.25, holdout= ''):
       #print('test@@@@@@', len(feature_sets), type(feature_sets))
        (label,tax,pool,pool_tax,baltimore) = subset
        self.label = label
        self.tax = tax
        self.pool = pool
        self.pool_tax = pool_tax
        self.baltimore = baltimore
        self.holdout = holdout
        if self.holdout:
            self.holdout = holdout
            self.ds = self.ho_get_class_lists(subset,vh)
            print ('holdout ', holdout, len(self.ds))
        else:   
            self.ds = self.get_class_lists(subset,vh)
        self.removeMissing()
        self.tt_split=tt_split
        if self.tt_split >0:
            self.trn_tst_split()
        print(self.ds.head())
        self.fs = feature_sets
        if feature_sets:
            print('Adding feature sets ',feature_sets)
            self.addFeatureSets(feature_sets)
        else:
            self.fs =[]
            
        
    def printX(self):
        print('this is an instance method')
        
    def get_class_lists(self,ss,vhDB):
        """ returns a dataframe virus,y,refseqs"""    
    
        (label,label_tax,pool,pool_tax,baltimore) = ss
        if baltimore =='all':
            baltimore = 'NA'  # will pick up all classes (except satallites)
        viruses = vhDB.viruses
        hosts = vhDB.hosts
        refs =[]
        # Get a list of all the viruses in the labelled class and one for the rest of the pool
        bc = [[],[]] 
        for v, vd in viruses.items() :
            if baltimore in vd['class']:
                host_labels = [hosts[h][label_tax] for h in vd['hosts'] if hosts[h][pool_tax] == pool] 
                if len(host_labels) > 0:
                # check for missing references or duplicates( phage in VHBB twice with same host at different tax levels)
                    if[r  for r in vd['refseqs']if r not in refs]: # ref not in list already
                        if label in host_labels:
                            bc[0].append ((v,1))
        
        # Add positive vireuses first then only add viruses with new refsseqs
        for v, vd in viruses.items() :
            if baltimore in vd['class']:
                host_labels = [hosts[h][label_tax] for h in vd['hosts'] if hosts[h][pool_tax] == pool] 
                if len(host_labels) > 0:
                # check for missing references or duplicates( phage in VHBB twice with same host at different tax levels)
                    if[r  for r in vd['refseqs']if r not in refs]: # ref not in list already
                        if label not in host_labels:
                            bc[1].append((v,0)) 

        # Get a random sample for each class of size n , the size of smallest class                    
        n = min(len(bc[0]),len(bc[1]),200)
        #n = min(len(bc[0]),len(bc[1]),10) #  for testing
        datas =  []
        for clss in bc: 
            datas.extend(random.sample(clss, n))
        
        dataset ={}
        for v,y in datas:
            dataset [v]= {'y':y, 'refseqs':viruses[v]['refseqs']} 
        df = pd.DataFrame(dataset).T
        df.index.name = 'virus'
        df.reset_index(inplace=True)
        return df
    
    
    def remove_close_viruses (self,ho, not_ho,ani,viruses):
        ho_refs = {r:v  for  v in ho for r in viruses[v]['refseqs']if r in ani.index}
        not_ho_refs = {r:v  for v in not_ho for r in viruses[v]['refseqs']if r in ani.index}
        too_close = []
        for ref in not_ho_refs:
            if(ani.loc[ref][ho_refs].values > 75).any():
                too_close.append(not_ho_refs[ref] ) 
        not_ho2 = [ v for  v in not_ho if v not in too_close ]
        return not_ho2

    def check_4_close_viruses(self,ho_pos, ho_neg,not_ho_pos,not_ho_neg,viruses):
        numbers = [len(not_ho_pos), len( not_ho_neg)]
        print('before ani -not_ho pos',len(not_ho_pos), 'neg',len(not_ho_neg)  )
        ani = pd.read_csv('ani.csv',index_col=0)
        not_ho_pos = self.remove_close_viruses (ho_pos, not_ho_pos,ani,viruses)
        not_ho_neg = self.remove_close_viruses (ho_neg, not_ho_neg,ani,viruses)
        print('after ani -not_ho pos',len(not_ho_pos), 'neg', len(not_ho_neg) )
        numbers.extend([len(not_ho_pos),len(not_ho_neg)])
        print (f'{self.label}_{self.holdout}',numbers)
#         csvfile = 'ANI_numbers.csv'       
#         with open(csvfile, 'a') as csvfile: 
#             writer = csv.writer(csvfile, delimiter=',',)
#             writer.writerow([f'{self.label}_{self.holdout}',numbers]) 
        
        return  not_ho_pos,not_ho_neg 

    def create_df (self,ho_pos, ho_neg,not_ho_pos,not_ho_neg,viruses):
        dataset ={}
        print('ho_pos',len(ho_pos),'ho_neg',len(ho_neg),'not_ho_pos',len(not_ho_pos),'not_ho_neg',len(not_ho_neg))
        for v in ho_pos:
            
            dataset[v]= {'y':1,'refseqs':viruses[v]['refseqs'], 'trn/tst': 'test'}
        for v in ho_neg:
            dataset[v]= {'y':0,'refseqs':viruses[v]['refseqs'], 'trn/tst': 'test'}
        for v in not_ho_pos:
            dataset[v]= {'y':1,'refseqs':viruses[v]['refseqs'], 'trn/tst': 'train'}
        for v in not_ho_neg:
            dataset[v]= {'y':0,'refseqs':viruses[v]['refseqs'], 'trn/tst': 'train'}

        df = pd.DataFrame(dataset).T
        df.index.name = 'virus'
        df.reset_index(inplace=True)

        return df

#     def get_class_lists(self,subset,vhDB):
#""" Attempting to use ANI aware method but need to decide when to remove close viruses
#"""
#         print('in new get class lists')
#         (label,label_tax,pool,pool_tax,balt) = subset
#         viruses = vhDB.viruses
#         hosts = vhDB.hosts
        
#         data_lists = {'training':[[],[]],'test':[[],[]]}
#         # Get a list of all the viruses in the labelled class and the rest of the pool
#         baltimore = 'NA'  # all slasses (except satallites)
#         tst_pos =[]
#         tst_neg =[]
#         trn_pos = []
#         trn_neg = []
#         all_data = []
#         for v, vd in viruses.items() :
#             if  baltimore in vd['class']: 
#                 host_labels = [hosts[h][label_tax] for h in vd['hosts'] if hosts[h][pool_tax] == pool]
#                 if host_labels:
#                     if label in host_labels:
#                         trn_pos.append(v)
#                     else:
#                         trn_neg.append(v)
# if[r  for r in vd['refseqs']if r not in refs]: # ref not in list already



#         max_tst = min (len(tst_pos),len(tst_neg),50)
#         tst_pos =    (random.sample(tst_pos, max_tst)) 
#         tst_neg =    (random.sample(tst_neg, max_tst))

#         # Remove viruses from training that are > 75% ANI to holdout viruses
#         trn_pos,trn_neg = self.check_4_close_viruses(tst_pos, tst_neg,trn_pos,trn_neg,viruses)

#         max_trn = min (len(trn_pos2),len(trn_neg),200)
#         trn_pos2 =    (random.sample(trn_pos, max_trn)) 
#         trn_neg2 =    (random.sample(trn-_neg, max_trn))

#         df = self.create_df (tst_pos, tst_neg,trn_pos,trn_neg,viruses)
#         return df
    
    def ho_get_class_lists(self,subset,vhDB):

        (label,label_tax,pool,pool_tax,balt) = subset
        viruses = vhDB.viruses
        hosts = vhDB.hosts
        
        data_lists = {'training':[[],[]],'test':[[],[]]}
        # Get a list of all the viruses in the labelled class and the rest of the pool
        baltimore = 'NA'  # all slasses (except satallites)
        ho_pos =[]
        ho_neg =[]
        not_ho_pos = []
        not_ho_neg = []
        all_data = []
        for v, vd in viruses.items() :
            if  baltimore in vd['class']: 
                host_labels = [hosts[h][label_tax] for h in vd['hosts'] if hosts[h][pool_tax] == pool]
                if host_labels:
                    if  vd['family']!= self.holdout:
                        ds_name ='not_ho'
                        if label in host_labels:
                            not_ho_pos.append(v)
                        else:
                            not_ho_neg.append(v)

                    else:
                        ds_name = 'ho'
                        if label in host_labels:
                            ho_pos.append(v)
                        else:
                            ho_neg.append(v)


        max_ho = min (len(ho_pos),len(ho_neg),20)
        ho_pos2 =    (random.sample(ho_pos, max_ho)) 
        ho_neg2 =    (random.sample(ho_neg, max_ho))

        # Remove viruses from training that are > 75% ANI to holdout viruses
        not_ho_pos2,not_ho_neg2 = self.check_4_close_viruses(ho_pos2, ho_neg2,not_ho_pos,not_ho_neg,viruses)

        max_not_ho = min (len(not_ho_pos2),len(not_ho_neg2),200)
        not_ho_pos2 =    (random.sample(not_ho_pos2, max_not_ho)) 
        not_ho_neg2 =    (random.sample(not_ho_neg2, max_not_ho))

        df = self.create_df (ho_pos2, ho_neg2,not_ho_pos2,not_ho_neg2,viruses)
        return df
    
    def removeMissing(self):
        missing = []
        for refs in self.ds['refseqs']:
            for ref in refs:
                for  fp,ext in self.filepaths:
                    f = fp/ (ref.strip() +ext )
                    if (not f.exists()) or (not f.read_text()): #if not filename.exists():
                        missing.append(ref)
        
        # remove viruses with at least 1 missing files
        refs = {r:i  for i, rs in enumerate (self.ds ['refseqs'])for r in rs }
        todrop = set([refs[r]for r in missing])
        self.ds.drop(todrop,axis=0 ,inplace = True)
        return 
         
                    
    def trn_tst_split(self,  validation = False):
        """ use sk learn to get stratified train test split and add a col for trn/tst/val
            need to deal with validation set 
        """
        
        train,test = train_test_split(self.ds, test_size=self.tt_split, random_state=0, stratify=self.ds[['y']]) #sk_learn
        print ('train test split', len(train),len(test))
        train['trn/tst'] = train.apply(lambda row: 'train', axis = 1)
        test['trn/tst'] = test.apply(lambda row: 'test', axis = 1)
        #val['trn/tst'] = datasets.apply(lambda row: 'val', axis = 1)
        df = pd.concat([train,test],axis=0,sort=False)
        #pd.merge(df,df2[['Key_Column','Target_Column']],on='Key_Column', how='left')
        self.ds = pd.merge (self.ds, df[['virus' ,'trn/tst']], on='virus')
        return 
        
    def addFeatureSets(self,f_sets):
        for fs in f_sets:
            print('adding fs',fs)
            self.fs.append(FeatureSet(self.ds, fs,self.filepaths))
            
    def results2CSV(self,results, subset, csvfile):
        (label,label_tax,pool,pool_tax,balt) = subset
        #results = { k : round(v, 3) for k,v in results.items()}
        results['label'] = self.label
        results['label_tax']= self.tax
        results['pool_label']= self.pool
        results['pool_tax']= self.pool_tax
        results['Baltimore'] = self.baltimore
        results['holdout'] = self.holdout
        results['N'] = len(self.ds)
        print(results)

        if os.path.isfile(csvfile):
           
            with open(csvfile, 'a') as csvfile:
                fieldnames = ['label','label_tax','pool_label','pool_tax','holdout','Baltimore', 'N' , 'features','k','AUC',\
                             'Acc', 'Spec', 'Sens', 'Prec', 'TP', 'TN', 'FP', 'FN']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writerow(results) 
        else:
            with open(csvfile, 'a') as csvfile:
                print ( 'new file',csvfile)
                fieldnames = ['label','label_tax','pool_label','pool_tax','holdout','Baltimore', 'N' , 'features','k','AUC',\
                             'Acc', 'Spec', 'Sens', 'Prec', 'TP', 'TN', 'FP', 'FN']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(results)  
