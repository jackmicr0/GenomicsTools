# GenomicsTools

A selection of tools for genomics data

- Retrieve_with_EaZeEeeee.py allows you to retrieve sequences and metadata as if you were searching on the ncbi search bar. outputs a fasta of all returned sequences and a csv file with metadata associated. 

- blastn_with_readable_descriptions.py allows local blasting but will add useful descriptions to output csv

- Primer_scheme_gap_finder_mean.py requires bed file for primer scheme and output of samtools depth as input. Will get mean depth across primer pair and amplicon.

- Primer_scheme_gap_finder_minimum_depth.py requires bed file for primer scheme and output of samtools depth as input. Will get minimum depth across primer pair and amplicon.
    - most useful if depth files are made with seperate raw data for each primer scheme pool as this prevents overlapping reactions across amplicons.
    - additional options for visualisation of depths
