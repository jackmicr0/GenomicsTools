#!/usr/bin/python3


## Jack Crook

import sys, os
import pandas as pd
import subprocess
import argparse
import argcomplete
from Bio import SeqIO, Entrez


## user args
argParser = argparse.ArgumentParser(description="Python tools to get blastn matches for fasta file input, returns NCBI description matching accession of target matches. Requires local blasting db e.g from makeblastdb tool.")
argParser.add_argument("-f", "--fasta", type=str, help="Fasta file of sequences of interest")
argParser.add_argument("-db", "--database", type=str, help="path to local database (incl prefix of database files)")
argParser.add_argument("-o", "--output_csv", type=str, help="output file name with extension")
argParser.add_argument("-e", "--email", type=str, help="email for ncbi tools")

argcomplete.autocomplete(argParser)

args = argParser.parse_args()

### function to get descriptions from accessions :)

def retrieve_record(accession): # Function to retrieve records from NCBI using accession numbers
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")         # Fetch the record from NCBI by accession number
        record = SeqIO.read(handle, "fasta")
        description = record.description
        handle.close()
        return description
    except Exception as e:
        print(f"Error retrieving accession {accession}: {str(e)}")
        return None


# Define your input FASTA file and output CSV file names
#input_fasta = sys.argv[1]
#db_prefix =sys.argv[2] ### db_prefix needs to be complete path to blast db - incl prefix for db files
#output_csv = sys.argv[3] 
#email = sys.argv[4]

Entrez.email=args.email

# Define the BLAST command with the appropriate database and options
blast_command = [
    "blastn",                  # BLASTn command
    "-db", args.database,    # Replace 'your_database' with your actual database name or path
    "-query", args.fasta,     # Input FASTA file
    "-outfmt", "6",           # Output format (CSV)
    "-out", args.output_csv, # Output CSV file
    "-max_target_seqs", "5" ## 5 top blast hits 
]

# Run the BLASTn command using subprocess
try:
    subprocess.run(blast_command, check=True)
    print("BLASTn search completed successfully.")

except subprocess.CalledProcessError as e:
    print("Error running BLASTn:", e)

    exit(1)


df = pd.read_csv(args.output_csv,sep="\t",names=["query id",
        "target id",
        "percentage id",
        "aln length",
        "mismatch",
        "gap openings",
        "query start",
        "query end",
        "target start",
        "target end",
        "evalue (significance)",
        "bitscore"])



print("Just getting human readable time saving nuggets of gold..")

df["Record"] = df["target id"].apply(retrieve_record)

df=df[["query id",
        "target id",
	"Record",
        "percentage id",
        "aln length",
        "mismatch",
        "gap openings",
        "query start",
        "query end",
        "target start",
        "target end",
        "evalue (significance)",
        "bitscore"]]


df.to_csv(args.output_csv,sep=",",index=False)

print("Commence virus hunting! woo.. ")
