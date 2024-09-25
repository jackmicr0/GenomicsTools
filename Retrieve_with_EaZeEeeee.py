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

### get accessions, sequences and metadata for your perusal.. 

argParser = argparse.ArgumentParser(description="Retrieve metadata and sequences with EAZEEeeeeeEEEEE - Use search query as you would if on NCBI directly and get sequences and metadata for that search")
argParser.add_argument("-q","--query_search", type=str, help="Query search for finding genomes e.g. SARS AND CoV NOT bacteria NOT unverified", nargs="+")
argParser.add_argument("-of","--output_fasta", type=str,help="Fasta file name with extension, e.g. sequences.fasta")
argParser.add_argument("-oc","--output_csv",type=str,help="Csv output file name with extension e.g. output.csv")
argParser.add_argument("-e","--email", type=str,help="Email address for EFetch tools")

args = argParser.parse_args() ### define args

## handle query terms 

query_words = args.query_search
query = " ".join(query_words)

### functions for running script 

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
    
def fetch_sequences(accessions):
    records = []
    for accession in accessions:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta")
        record = SeqIO.read(handle, "fasta")
        records.append(record)
        handle.close()
    return records

def create_sequence_dataframe(records):
    df = pd.DataFrame({
        'Accession': [record.id for record in records],
        'Description': [record.description for record in records],
        'Sequence': [str(record.seq) for record in records]
    })
    df['Corrected_accession'] = df['Accession'].str.split('.').str[0]
    return df[['Accession', 'Description', 'Sequence', 'Corrected_accession']]

def fetch_cds_and_metadata(accessions):
    cds_list = []
    feature_list = []
    handle = Entrez.efetch(db="nucleotide", id=accessions, rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "genbank")
    
    for record in records:
        accession = record.name
        for feature in record.features:
            if feature.type == "CDS":
                cds_list.append({
                    'Accession': accession,
                    'start': feature.location.start,
                    'end': feature.location.end,
                    'strand': feature.location.strand
                })
            if feature.type == "source":
                feature_info = {
                    'Accession': accession,
                    'isolation_source': feature.qualifiers.get('isolation_source', [np.nan])[0],
                    'isolate': feature.qualifiers.get('isolate', [np.nan])[0],
                    'geo_loc_name': feature.qualifiers.get('geo_loc_name', [np.nan])[0],
                    'host': feature.qualifiers.get('host', [np.nan])[0],
                    'segment': feature.qualifiers.get('segment', [np.nan])[0],
                    'mol_type': feature.qualifiers.get('mol_type', [np.nan])[0],
                    'collection_date': feature.qualifiers.get('collection_date', [np.nan])[0]
                }
                feature_list.append(feature_info)
    
    handle.close()
    
    df_cds = pd.DataFrame(cds_list)
    df_features = pd.DataFrame(feature_list)
    
    return df_cds, df_features

def clean_dataframe(df):
    # Apply to each column, ensuring it applies to string columns only
    df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))
    df = df.apply(lambda col: col.map(lambda x: re.sub(r'[-/()|.:;""\' ]', '_', x) if isinstance(x, str) else x))
    df = df.apply(lambda col: col.map(lambda x: x.replace("__", "_") if isinstance(x, str) else x))
    df = df.apply(lambda col: col.map(lambda x: x.strip('_') if isinstance(x, str) else x))
    return df

def fill_missing_data(df):
    df.fillna({
        "isolation_source": 'NoIsolationSource',
        "collection_date": 'NoCollectionDate',
        "isolate": 'NoIsolate',
        "geo_loc_name": 'NoCountry',
        "host": 'NoHost',
        "mol_type": 'NoMolType',
        "start": 'NoStartCoordinate',
        "end": 'NoEndCoordinate',
        "segment": 'NoSegmentOrNotSegmented'
    }, inplace=True)
    return df

def create_fasta_headers(df):
    df["FastaID"] = df.index + "|" + df["isolate"] + "|" + df["host"] + "|" + df["isolation_source"] + "|" + df["geo_loc_name"] + "|" + df["segment"] + "|" + df["collection_date"]
    return df

def process_accessions(accessions, output_fasta, output_csv):
    # Fetch sequences and metadata
    sequence_records = fetch_sequences(accessions)
    df_unfiltered = create_sequence_dataframe(sequence_records)
    
    # Fetch CDS and metadata
    df_cds, df_features = fetch_cds_and_metadata(accessions)
    
    # Merge DataFrames
    df_unfiltered.set_index("Corrected_accession", inplace=True)
    df_cds.set_index("Accession", inplace=True)
    df_features.set_index("Accession", inplace=True)
    
    merged_df = pd.merge(df_unfiltered, df_cds, left_index=True, right_index=True)
    merged_df = pd.merge(merged_df, df_features, left_index=True, right_index=True)
    
    # Clean and handle missing data
    merged_df = clean_dataframe(merged_df)
    merged_df = fill_missing_data(merged_df)
    
    # Create FASTA headers
    merged_df = create_fasta_headers(merged_df)
    
    # Write FASTA file
    with open(output_fasta, 'w') as fasta_file:
        for _, row in merged_df.iterrows():
            fasta_file.write(f">{row['FastaID']}\n{row['Sequence']}\n")
    
    # Write CSV file
    merged_df.to_csv(output_csv, sep=",")


if __name__ == "__main__":
    # Retrieve accessions using the provided search query
    accessions_out_noduplicates = get_accessions(query)

    if not accessions_out_noduplicates:
        print("No accessions found for the given search query.")
    else:
        process_accessions(accessions_out_noduplicates, args.output_fasta, args.output_csv)
        print(f"Results written to {args.output_fasta} and {args.output_csv}")
