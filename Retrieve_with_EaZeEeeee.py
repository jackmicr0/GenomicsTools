#!/usr/bin/env python
#
#Created on Wed Aug 16 11:40:31 2023

#@author: Jack.Crook

from Bio import Entrez
from Bio import SeqIO
import numpy as np
import pandas as pd
import re
import argparse

### get lovely accessions, sequences and metadata for your perusal.. 

argParser = argparse.ArgumentParser(description="Retrieve metadata and sequences with EAZEEeeeeeEEEEE - Use search query as you would if on NCBI directly and get sequences and metadata for that search")
argParser.add_argument("-q","--query_search", type=str, help="Query search for finding genomes e.g. SARS AND CoV NOT bacteria NOT unverified", nargs="+")
argParser.add_argument("-of","--output_fasta", type=str,help="Fasta file name with extension, e.g. sequences.fasta")
argParser.add_argument("-oc","--output_csv",type=str,help="Csv output file name with extension e.g. output.csv")
argParser.add_argument("-e","--email", type=str,help="Email address for EFetch tools")

args = argParser.parse_args() ### define args

## handle query terms 

query_words = args.query_search
query = " ".join(query_words)

### functions

def get_accessions(search_query):
    
    Entrez.email = args.email # Provide your email address

    handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()

    if "IdList" in record:
        id_list = record["IdList"]
        fetch_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="acc", retmode="text")
        acc_ids = fetch_handle.read().strip().split()
        fetch_handle.close()
        return acc_ids
    else:
        return []

search_query = query

accessions_out = get_accessions(search_query)
accessions_out_noduplicates=list(set(accessions_out)) ## no duplicates 
accessions_for_retrieval = accessions_out_noduplicates


### get cds start and end and metadata for each accession 

records = [] ## empty list for saving records 
df=pd.DataFrame() ## start empty df 
for accession in accessions_out_noduplicates: ## loop accessions - get fasta 
    handle = Entrez.efetch(db="nucleotide",\
                           id=accession,\
                           rettype="fasta")
    
    record = SeqIO.read(handle, "fasta") # get sequence 
    
    records.append(record)

    
df['Accession'] = [record.id for record in records]
df['Description'] = [record.description for record in records]
df['Sequence'] = [str(record.seq) for record in records]
df_unfiltered=df[['Accession','Description','Sequence']]
df_unfiltered['Corrected_accession']=df_unfiltered['Accession'].str.split(".").str[0] 

################## METADATA

cds_list=[] ## initiate end list 


handle = Entrez.efetch(db="nucleotide",\
                      id=accessions_for_retrieval,\
                      rettype="gb",\
                      retmode="text")

records = SeqIO.parse(handle, "genbank")

for record in records:
   
    for feature in record.features:
        
        if feature.type == "CDS":
            cds = {} 
            cds["Accession"] = record.name         
            cds["start"] = feature.location.start                 
            cds["end"] = feature.location.end
            cds["strand"] = feature.location.strand
                
            cds_list.append(cds)   
                
handle.close()

df_cds = pd.DataFrame(cds_list) ## df with start and end coordinates 

feature_list=[]
handle = Entrez.efetch(db="nucleotide",\
                      id=accessions_for_retrieval,\
                      rettype="gb",\
                      retmode="text")
records = SeqIO.parse(handle, "genbank")

for record in records:   
    for feature in record.features:
        
        if feature.type == "source":                 
                feature_info = {}                 
                feature_info["Accession"] = record.name                 
                if "isolation_source" in feature.qualifiers: 
                    feature_info["isolation_source"] = \
                    feature.qualifiers["isolation_source"][0]                 
                else:                     
                    feature_info["isolation_source"] = float('nan')
                
                if "isolate" in feature.qualifiers:                     
                    feature_info["isolate"] = \
                    feature.qualifiers["isolate"][0]
                else:                     
                    feature_info["isolate"] = float('nan')                 
                if "country" in feature.qualifiers:                     
                    feature_info["country"] = \
                    feature.qualifiers["country"][0]
                else:                  
                    feature_info["country"] = float('nan') 
                if "host" in feature.qualifiers:       
                        feature_info["host"] = \
                        feature.qualifiers["host"][0]
                else:                     
                    feature_info["host"] = float('nan')   
                if "segment" in feature.qualifiers: 
                    feature_info["segment"] = \
                    feature.qualifiers["segment"][0]   
                else:                   
                    feature_info["segment"] = float('nan') 
                if "mol_type" in feature.qualifiers:       
                        feature_info["mol_type"] = \
                        feature.qualifiers["mol_type"][0]
                else:                     
                    feature_info["mol_type"] = float('nan')
                if "collection_date" in feature.qualifiers:       
                        feature_info["collection_date"] = \
                        feature.qualifiers["collection_date"][0]
                else:                     
                    feature_info["host"] = float('nan')
                 
                
                feature_list.append(feature_info)    
                
handle.close()

df_features = pd.DataFrame(feature_list) 

df_unfiltered.set_index("Corrected_accession",inplace=True)

df_cds.set_index("Accession", inplace=True)
df_features.set_index("Accession", inplace=True)

merge = pd.merge(df_unfiltered, df_cds, right_index=True, left_index=True).merge(df_features, right_index=True, left_index=True)

merge = merge.applymap(lambda x: x.strip() if isinstance(x, str) else x)
merge = merge.applymap(lambda x: re.sub(r'[-/()|.:;""'' ]', '_', x) if isinstance(x, str) else x)
merge = merge.applymap(lambda x: x.replace("__", "_") if isinstance(x, str) else x)
merge = merge.applymap(lambda x: x.replace('""', "") if isinstance(x, str) else x)
merge = merge.applymap(lambda x: x.rstrip("_") if isinstance(x, str) else x)
merge = merge.applymap(lambda x: x.lstrip("_") if isinstance(x, str) else x)

merge["isolation_source"] = merge["isolation_source"].replace(np.NaN, 'NoIsolationSource')
merge["collection_date"] = merge["collection_date"].replace(np.NaN, 'NoCollectionDate')
merge["isolate"] = merge["isolate"].replace(np.NaN, 'NoIsolate')
merge["country"] = merge["country"].replace(np.NaN, 'NoCountry')
merge["host"] = merge["host"].replace(np.NaN, 'NoHost')
merge["mol_type"] = merge["mol_type"].replace(np.NaN, 'NoMolType')
merge["start"] = merge["start"].replace(np.NaN, 'NoStartCoordinate')
merge["end"] = merge["end"].replace(np.NaN, 'NoEndCoordinate')
merge["segment"] = merge["segment"].replace(np.NaN, 'NoSegmentOrNotSegmented')

merge["FastaID"]=merge.index + "|" + merge["isolate"] + "|" + merge[
    "host"] + "|" + merge[
    "isolation_source"] + "|" + merge[
    "country"] + "|" + merge[
    "segment"] + "|" + merge[
    "collection_date"
]
        
## make fasta 

output_file = args.output_fasta

# Open the file in write mode and iterate over DataFrame
with open(output_file, 'w') as fasta_file:
    
    for index, row in merge.iterrows():
        header = '>' + row['FastaID']  # Create FASTA header
        sequence = row['Sequence']    # Extract sequence

        # Write header and sequence to the file
        
        fasta_file.write(header + '\n')
        fasta_file.write(sequence + '\n')

        
merge.to_csv(args.output_csv, sep=",")
