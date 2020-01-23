import numpy as np
from Bio import SeqIO
import sys
from itertools import product
from collections import defaultdict, Counter
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.preprocessing import StandardScaler

class FeatureSet(object):
    
    genome_reps = ['DNA','AA','PC','Domains']
    pc_bins = {'A':'t','G':'t','V':'t','I':'u','L':'u','F':'u','P':'u',
                'M':'v','S':'v','T':'v','Y':'v', 'H':'w','N':'w','Q':'w','W':'w',
                'R':'x','K':'x','D':'y','E':'y','C':'z'
              }
  
    na_dict = 'acgt'
    aa_dict = 'ARNDCEQGHILKMFPSTWYV'
    pc_dict = 'tuvwxyz' # names of physiochem. bins
    symbol_dicts = [na_dict,aa_dict,pc_dict, {}]
   
    def __init__(self,ds,f_s,filepaths):  
        self.name    = f_s
        self.feature = f_s.split('_')[0]
        self. k      = int(f_s.split('_')[1])
        self.index   =  self.genome_reps.index (self.feature)
        self.symbols = self.symbol_dicts[self.index]
        print('FS get sequences for length ds', len(ds))
        seqs = self.get_sequences(ds,filepaths)
        print('FS get feature names', f_s, len(seqs))
        self.feature_names = self.get_feature_names(seqs)
        print('FS get X','len fs',len(self.feature_names))
        X,X_trn, X_tst,K_trn,K_tst = self.get_X (ds,seqs)
        self.X = X
        self.X_trn = X_trn
        self.X_tst = X_tst
        self.K_trn = K_trn
        self.K_tst = K_tst
   
    def __str__(self):
        return f'Feature Set {self.name}'
    
    def get_sequences(self, df,filepaths):
        fp,ext = filepaths[self.index]
        sequences = defaultdict(str)
        missing_files = []
           
        for i, row in df.iterrows():
            v=(row['virus'])
            seqs = []
            for ref in row['refseqs']:
                filename = (fp/ref.strip()).with_suffix(ext)
                try:
                    with open(filename, 'rt') as handle:
                        if self.feature != 'Domains':
                            for record in SeqIO.parse(handle, "fasta"):                   
                                  seqs += str(record.seq)
                        else:
                            line = handle.read()
                            seqs.append (','.join (line.split())) 
                except IOError:
                     missing_files.append(v)
            
            if self.feature == 'PC': # bin PC
                pc_seqs = ''.join([self.pc_bins.get(s,'X') for s in seqs ])
                sequences [v] = ''.join(pc_seqs)   
            elif  self.feature == 'Domains':   
                sequences [v] = ','.join(seqs)
            else:
                sequences [v] = ''.join(seqs)
#             if i%100  == 0:
#                 print(i, filename,v,len(sequences[v]),sequences[v][:40])
        #print(missing_files)
        return sequences
    
    def get_X(self,ds,seqs  ):

        X = self.extract_kmers(ds,seqs)
        mask = ds['trn/tst']=='train'
        X_train = X[mask,:]
        X_test = X[~mask,:]
        

        scaler = StandardScaler()
        X_trn_std = scaler.fit_transform(X_train)
        X_tst_std = scaler.transform(X_test)
        gram_train = np.dot(X_trn_std,X_trn_std.T)
        gram_tst = np.dot(X_tst_std, X_trn_std.T)
        return X,X_trn_std, X_tst_std,gram_train,gram_tst
    
    def extract_kmers (self, df,sequences):
      
        seq_index = []
       
        X = np.zeros((len(df) ,len(self.feature_names)))
        print(len(df),len(sequences) ,len(self.feature_names), np.shape(X))
        for i,(ind, row) in enumerate(df.iterrows()):
            v=row['virus']
            seq = sequences[v]
            if seq:
                if self.feature == 'Domains':
                    for dom in seq.split(','):
                        j = self.feature_names[dom]
                        X[i,j] += 1 
                else:
                    X[len(seq_index)] = self.get_word_frequ (seq)
                seq_index.append(v)
        # normalise for length of genome(s)
        XN = np.divide(X,X.sum(axis = 1)[:,None])
        return XN#,seq_index,f_index

    def get_feature_names(self,sequences):
        if self.feature == 'Domains': 
            all_kmers = {dom for s in sequences.values() for dom in s.split(',')}
        else:
            all_kmers =[''.join(c) for c in product(self.symbols, repeat=self.k)]
        return {kmer:i for i,kmer in enumerate (all_kmers)}
    
    def get_word_frequ(self,sequence):
        """returns the word frequency array for all genomes for a virus """
        word_freqs = np.zeros(shape =len(self.feature_names),dtype=np.int32)
        for i in range(len(sequence)-self.k+1):
            word = sequence[i:i+self.k]
            if self.check_word(word):
                kmer_index = self.feature_names[word]
                word_freqs [ kmer_index] += 1
                # If this for DNA/RNA then add the  compliment
                if self.feature=='DNA':
                     word_freqs [self.feature_names[self.compliment(word)]] +=1   
        return word_freqs
   
    def check_word(self,mystring):
        """Check for illegal characters """
        return all(c in self.symbols for c in mystring)

    def compliment(self,dna):
        """Return the reverse compliment of a sequence """
        bp = { 'a':'t', 't':'a', 'c':'g','g':'c'}
        revdna =''    
        for i in range(len(dna)-1,-1,-1):
             revdna += bp[dna[i]]
        return (revdna)        