import pandas as pd
import numpy as np
import re
import pickle
from collections import defaultdict
from ete3 import NCBITaxa



class VHDB:

    """ Class to hold all the information from the Virus Host Database """



    def __init__(self, inputfile):

        self.viruses = {}
        self.hosts = defaultdict(dict)
        self.subspecies =defaultdict(dict)
        self.refseq2taxid = {}
   
        self.loadVHDB(inputfile)
        #self.ncbi = NCBITaxa()
        
    def loadVHDB(self,inputfile):

        no_lineage = []
        lines = 0
        missing_host = 0
        
        vh_df = pd.read_csv(inputfile,na_filter=False)
        ncbi = NCBITaxa()

        for index, row in vh_df.iterrows():  
            lines += 1
            virus = str(row["virus tax id"])
            host = str (row["host tax id"])
            lineage = str(row["virus lineage"])
            

    # for all the entries that are Viruses(not Viroids) and have a host

            #if host == '' or host == '1' or  re.search (r'Viroid',lineage):
            if host == '' or host == '1' or re.search (r'Viroid',lineage): 
                missing_host  += 1
            else:
        
                self.hosts[host] = self.get_desired_ranks(host,ncbi)

                if not virus  in self.viruses:         
                    self.viruses[virus] =  self.get_metadata(row,ncbi) 

                else:  
                # virus already in dictionary add new host and new refs 
                    if not host in self.viruses[virus]['hosts']:                      
                        self.viruses[virus]['hosts'].append(host)
                    for r in str(row["refseq id"]).split(','):
                        if r not in self.viruses[virus]['refseqs']:
                            self.viruses[virus]['refseqs'].append(r) 


        print (f' strain level viruses {len(self.subspecies)},  viruses with no lineage info {len(no_lineage)}')
        print (f' number of viruses  {len(self.viruses)} from {lines} lines of csv')
        print (f' number of hosts {len (self.hosts)} and {missing_host} missing hosts')
            
                
        return 


    
    def get_desired_ranks(self,taxid, ncbi ):
        desired_ranks = ['superkingdom','kingdom', 'phylum', 'class', 'order', 'family', 'genus','species' ]
        try:
            lineage = ncbi.get_lineage(taxid)
            tax2ranks = ncbi.get_rank(lineage)
            taxid2name = ncbi.get_taxid_translator(tax2ranks)
            ranks2tax = {rank:taxid for taxid, rank in tax2ranks.items() if rank in desired_ranks}
            ranks2name = {rank:taxid2name[txid] for txid, rank in tax2ranks.items()}

            rank =  tax2ranks[int(taxid)]
            if rank in desired_ranks:
                pass
            elif  'sub' in rank:      
                taxid = ranks2tax[rank[3:]]
            else:  # rank not in desired_ranks
                #print ('rank not in desired_ranks',taxid,rank)
                if len (lineage )>1:
                    i = 1
                    while rank not in desired_ranks  :
                        taxid = lineage[-i]
                        rank =  tax2ranks[int(taxid)]
                        i += 1
                    taxid = ranks2tax[rank]

                else:
                    # host taxid = 1 
                    pass
                    print ('host',taxid,rank, taxid2name)
                    
        except:
             print(f'issue with taxid,{taxid}')
       
        tax_dict = {rank: ranks2name.get(rank, 'not assigned') for rank in desired_ranks}
        
        # for Bacteria and Archeae
        if tax_dict['kingdom']== 'not assigned':
            tax_dict['kingdom'] = tax_dict['superkingdom']
           
        return tax_dict
 
    def get_virus_ranks(self,taxid, ncbi ):
        
        desired_ranks = ['superkingdom','phylum', 'class', 'order', 'family', 'genus','species' ]
        lineage = ncbi.get_lineage(taxid)
        tax2ranks = ncbi.get_rank(lineage)
        rank2tax = {rank:tax for tax,rank in tax2ranks.items()}
        taxid2name = ncbi.get_taxid_translator(tax2ranks)
        ranks2name = {rank:taxid2name[txid] for txid, rank in tax2ranks.items()}
#         print (f'taxid {type(taxid)} lineage{lineage} \n  tax2ranks {tax2ranks} \n \
#         taxid2name {taxid2name} \n ranks2tax  {rank2tax}')
        
        # Dealing with viruses added at strain level getting the rank of the virus taxid should be species
        
        row_dict ={rank: ranks2name.get(rank, 'not assigned') for rank in desired_ranks}
        try:
            #print (tax2ranks[int(taxid)])
            rank =  tax2ranks[int(taxid)]
            if rank == 'no rank' :
                row_dict.update({'ss':taxid2name[int(taxid)] , 'species_taxid':rank2tax['species']})
            elif rank != 'species':
                print(f'****AGHHHH  {taxid}, rank {rank}')
        except:
             print(f'issue with taxid,{taxid},{ranks2name}')
 
        return row_dict

    
    def get_metadata(self,row,ncbi):
        """ create dictionry entry for a virus from pandas row """
        
        taxid = str(row["virus tax id"])
        tax_dict = self.get_virus_ranks(taxid, ncbi )
        balt =  self.getBaltimore(row)
        host =  [str(row['host tax id'])]
        refseqs = str(row["refseq id"]).split(',')
        
       
        if 'ss' in tax_dict:
            species_taxid = tax_dict['species_taxid']
            self.subspecies[species_taxid][taxid] = {'hosts':host,'refseqs':refseqs}
            
        tax_dict.update({'class':balt,'refseqs': refseqs, 'hosts': host}) 

        return tax_dict
    
    def getBaltimore (self,row):
    
        baltimore = [ 'dsDNA', 'ssDNA','DNA', 'dsRNA', 'RNA', '(+)ssRNA','(-)ssRNA','DNA_Retro', 'RNA_Retro', 'Satellites']
        balt_class = 'unassigned'
        lineage = str(row["virus lineage"])
        
        DNA_class = re.findall(r'\S*NA\b',lineage)
        if DNA_class:
            DNA_class = DNA_class[0].strip()
            if DNA_class in baltimore:
                balt_class = DNA_class
                
            elif 'RNA' in DNA_class :
                if re.search (r'positive',lineage):
                    balt_class = '(+)ssRNA'
                elif re.search (r'negative',lineage):
                    balt_class = '(-)ssRNA'
                else:
                     balt_class = 'RNA'
            
        elif re.search(r'Ortervirales',lineage):
            balt_class = 'RNA_Retro'
        elif re.search(r'Caulimoviridae', lineage)  or re.search('Hepadnaviridae',lineage):
            balt_class = 'DNA_Retro'
        elif re.search(r'Satellites',lineage):
            balt_class = balt_class = 'Satellites'
    
            
        return (balt_class)
######################################################################################
 
    

    
    
    
  
  
   

    def vhi(self):       
        v_h = np.zeros(shape=(n_viruses,n_hosts))
        for virus in virus_index:
            v = virus_index[virus]
            hosts = virus_metadata[virus]['hosts']
            for host in hosts:
                h = host_index[host]    
                v_h [v,h] = 1
                count += 1
        print ('There are {}interactions between {} viruses and {} hosts.'.format(count, n_viruses, n_hosts))
        return vhi


    def saveV_H(self,fileout):
        with  open( fileout, "wb" ) as fp:
            pickle.dump( self, fp  )
    
    @classmethod
    def loadV_H( cls,filein):
        with open(filein, 'rb') as fp:
            obj = pickle.load( fp)
        return obj
     



