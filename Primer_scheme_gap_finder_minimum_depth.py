#!/usr/bin/env python

'''

Author: Jack Crook
Email: Jack.crook@ukhsa.gov.uk

Find gaps in primer scheme data. Input files are depth files outputted from mapping raw fastq to reference used to generate the
scheme. Also need bed file with information for primer start and end site. bed file format doesnt matter as long as fwd and rv primers for each amplicon are one after the other.

How to format input files for using..

- Must map raw fastq to reference used for making primer scheme

- Sample_id must be first part of file name - before a "_" e.g. sampleid_restoffilename.depth (as sample_id is used for plot legend)

- If looking at individual pools, - option --primer_pairs must be a list of numbers wanted. --primer_pairs must also be used if x axis tick labels are wanted.

'''


import os, sys
import argparse
import pandas as pd
import argcomplete
import matplotlib.pyplot as plt
import glob
import numpy as np
import math
import warnings
warnings.filterwarnings("ignore")


'''

############# for testing in notebook #########################

bed = "depth10_cinderford_1kb_scheme.primer.bed"
pattern = "*.depth"
outdir_multi = "multi"
overall = "Overall.csv"
flower = "Flower.csv"
rlower = "Rlower.csv"
alower = "Alower.csv"
allcsv = "primer_amplicon.csv"
axhline = 50
per_plot = 10
xticks = []
log_condition = ""

############# Figures params ##################################

'''

plt.rcParams.update({"font.size":10})

################## ADD ARGUMENTS ##############################

description = """

Find gaps in primer scheme data - returns minimum depth in each region



Input files are depth files outputted from mapping raw fastq to reference used to generate the scheme and bed file with information for primer start and end site. bed file format needed:

	column 2 and 3 need to be start and end coords for primer.



How to format input files for using:

	- Must map raw fastq to reference used for making primer scheme.
	- Sample_id must be the first part of the file name (before a "_") e.g., sampleid_restoffilename.depth (as sample_id is used for plot legend).
	- If looking at individual pools, the --primer_pairs option must be a list of numbers wanted.	--primer_pairs must also be used if x-axis tick labels are wanted.





"""

argParser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
argParser.add_argument("-b","--bed", type=str, required=True, help="""Path to primer scheme bed file with coordinates for primers - make sure contains only pool being investigated""")
argParser.add_argument("-p","--pattern", type=str, required=True, help="Pattern to loop e.g. file extension, by default will look for .depth")

argParser.add_argument("-o","--outdir_multi_sample_plots", type=str, default='Plots', help="Make directory for multisample plot - default it 'Plots'")
argParser.add_argument("-a", "--all_output", type=str, default='Primer_and_amplicon_lower.csv', help="csv file name for lower output for depth across primer and amplicon - default is 'Primer_and_amplicon_lower.csv'")
argParser.add_argument("-f", "--flower_output", type=str, default='Forward_primer_lower.csv', help="csv file name for fwd primer lower output - default is 'Forward_primer_lower.csv'")
argParser.add_argument("-r", "--rlower_output", type=str, default='Reverse_primer_lower.csv', help="csv file name for rv primer lower output - default is 'Reverse_primer_lower.csv'")
argParser.add_argument("-am", "--alower_output", type=str, default='Amplicon_lower.csv', help="csv file name for amplicon lower output - default is 'Amplicon_lower.csv'")
argParser.add_argument("-l", "--axhline_for_plot", type=int, default=50, help="For visualisation: will put line on y axis of depth plots")
argParser.add_argument("-np", "--per_plot", type=int, default=10, help="Number of samples per plot")
argParser.add_argument("-pr", "--primer_pairs", type=str, default="", help="list of numbers of primer pairs for the pool of choice e.g. if primer pairs 1,3,5 are wanted - list should be 1,3,5")
argParser.add_argument("-s", "--scale", type=str, default="", help="write 'log' if want log scale, else deafult it not log")

args = argParser.parse_args() ### define args

### user specified strings for ease of use..

bed = args.bed
pattern = args.pattern
outdir_multi = args.outdir_multi_sample_plots
allcsv = args.all_output
flower = args.flower_output
rlower = args.rlower_output
alower = args.alower_output
axhline = args.axhline_for_plot
per_plot = args.per_plot
xticks = args.primer_pairs
log_condition = args.scale

if xticks and ',' in xticks:
    xticks = [int(x) for x in xticks.split(',')]
else:
    xticks = []

if not os.path.exists(outdir_multi): ### Make new directories if not in os.path
    os.makedirs(outdir_multi)

argcomplete.autocomplete(argParser) ### allow autocompletion


######################################## DEFINE FUNCTIONS FOR SCRIPT ########################################

def split_df(df, coordinates): ## Will split depth df into chunks according to coordinates in the bed file
    sub_dfs = []
    for start, end in coordinates:
        sub_dfs.append(df[start:end])
    return sub_dfs

def calculate_lower(sub_dfs, column_name): ## Will calculate the lowest depth across each amplicon
    lower_list = []
    for sub_df in sub_dfs:
        lower = round(sub_df[column_name].min())
        lower_list.append(lower)
    return lower_list

def create_new_df(mean_dict, primer_number): ## creates a new df with lower for each primer
    new_df = pd.DataFrame(mean_dict)
    new_df["primer_number"] = primer_number
    new_df

    return new_df

def create_and_explode_dataframe(amplicon_number_list, lower_list, sample_list, file_list): ## makes df with mean depths and formats
    df = pd.DataFrame()
    df['Amplicon_number'] = amplicon_number_list
    df['LowerDepth'] = lower_list
    df['Sample_id'] = sample_list
    df['File'] = file_list

    df_exp = df.apply(lambda col: col.explode() if col.dtype == 'object' else col)

    return df_exp

def plot_and_save_groups(df, per_plot, axhline, outdir, plot_output_prefix, x_tick_labels, log):

    groups = df.groupby("Sample_id")
    num_plots = math.ceil(len(groups) / per_plot)

    for i in range(num_plots):
        fig, ax = plt.subplots(figsize=(10, 5))

        for j, (name, group) in enumerate(list(groups)[i * per_plot:(i + 1) * per_plot]):
            ax.scatter(group["Amplicon_number"], group["LowerDepth"], label=name)

            if log == "log":
                plt.yscale("log", base=10)

            plt.xlabel("Amplicon number")
            plt.ylabel("Minimum depth")
            plt.axhline(axhline, color="k")

            x_min = min(group["Amplicon_number"])
            x_max = max(group["Amplicon_number"])
            x_tick_locs = range(int(x_min), int(x_max) + 1)
            plt.xticks(x_tick_locs, x_tick_labels)

            max_legend_col = 4
            num_entries = len(ax.lines)
            num_columns = int(np.ceil(num_entries / max_legend_col))
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, ncol=num_columns, loc="center left", bbox_to_anchor=(0.0, -0.3))

            plt.tight_layout()

        plt.savefig(os.path.join(outdir, f"{plot_output_prefix}_{i + 1}.png"), format="png", dpi=300)
        plt.close()
        plt.clf()

########################################## LOAD IN BED FILE ###################################################

with open(bed, "r") as bed_file:

    print(f"Bed file loaded: {bed}")

    bed_df=pd.read_csv(bed_file,sep="\t",names=["ref","start","end","pair","pairnum","fwdRev"])
    num_of_primer_pairs = int(len(bed_df.index) / 2) + 1 ## get number of primers
    repeated_numbers = np.repeat(np.arange(1, num_of_primer_pairs), 2) ## make repeated list to cover fwd + rv
    bed_df["Primer"]=repeated_numbers ## add to bed

    if not xticks:  # Check if xticks is an empty string
        mask = bed_df.index
    else:
        mask = bed_df['Primer'].isin(xticks)

        mask = bed_df['Primer'].isin(xticks) ## keep only those wanted
        bed_df = bed_df[mask]


bed_df["start"]=bed_df["start"].astype(int)
bed_df["end"]=bed_df["end"].astype(int)
bed_df["range"]=bed_df["end"]-bed_df["start"] ## get length of primer

primerfwd=bed_df.iloc[::2] ## forward primer rows - needed as each primer pair has 2 rows (fwd and rv)
primerrv=bed_df.iloc[1::2] ## reverse primer rows

fwd_coordinates = list(zip(primerfwd["start"], primerfwd["end"]))
rv_coordinates = list(zip(primerrv["start"], primerrv["end"]))
amp_coordinates = list(zip(primerfwd["end"], primerrv["start"]))
primers_and_amp_coords = list(zip(primerfwd["start"],primerrv["end"]))

amplicons = len(fwd_coordinates)

########################################## LOAD IN DEPTH FILE #########################

column_name = "depth"

########### loop directory and plot onto one plot ###############

# create lists for storage

fwd_lower_list = []
rv_lower_list = []
amplicon_lower_list = []
primer_amplicon_lower_list = []

fwd_master_amplicon_number_list = []
rv_master_amplicon_number_list = []
amplicon_master_number_list = []
primer_amplicon_number_list = []

fwd_master_sample_list = []
rv_master_sample_list = []
amplicon_master_sample_list = []
primer_amplicon_sample_list = []

fwd_master_file_list = []
rv_master_file_list = []
amplicon_master_file_list = []
primer_amplicon_master_file_list = []

fig, ax = plt.subplots(figsize=(15, 8)) # create ax object for plotting

files = glob.glob(f"*{pattern}")

for i, file in enumerate(files):

    print(f"Depth file: {file}")
    amplicon_num_list = np.arange(1, amplicons+1, 1) ## amplicon number being analysed
    file_temp = [file] * amplicons     # populate file list with N number of file labels

    fwd_master_file_list.append(file_temp)
    rv_master_file_list.append(file_temp)
    amplicon_master_file_list.append(file_temp)
    primer_amplicon_master_file_list.append(file_temp)

    filename=file.split("_")[0]
    print(f"Sample ID: {filename} (if not correct, make sure information before first '_' in filename is sampld id)")     # populate sample list with N number of sample labels
    sample_id = [filename] * amplicons

    fwd_master_sample_list.append(sample_id)
    rv_master_sample_list.append(sample_id)
    amplicon_master_sample_list.append(sample_id)
    primer_amplicon_sample_list.append(sample_id)

    df=pd.read_csv(file,sep="\t",names=["ref","base","depth"])

    fwd_sub_dfs = split_df(df, fwd_coordinates)     # subset df by coordinates given - as wm function loops through
    rv_sub_dfs = split_df(df, rv_coordinates)
    amp_sub_dfs = split_df(df, amp_coordinates)
    primer_amp_sub_dfs = split_df(df, primers_and_amp_coords)

    fwd_lower = calculate_lower(fwd_sub_dfs, column_name) ## calculate lower for each
    rv_lower = calculate_lower(rv_sub_dfs, column_name)
    amplicon_lower = calculate_lower(amp_sub_dfs, column_name)
    primer_amplicon_lower = calculate_lower(primer_amp_sub_dfs, column_name)

    fwd_lower_list.append(fwd_lower) ## append into list
    rv_lower_list.append(rv_lower)
    amplicon_lower_list.append(amplicon_lower)
    primer_amplicon_lower_list.append(primer_amplicon_lower)

    fwd_master_amplicon_number_list.append(amplicon_num_list)
    rv_master_amplicon_number_list.append(amplicon_num_list)
    amplicon_master_number_list.append(amplicon_num_list)
    primer_amplicon_number_list.append(amplicon_num_list)


## get exploded dfs

lower_fwd_df_exp = create_and_explode_dataframe(fwd_master_amplicon_number_list,fwd_lower_list,
                                               fwd_master_sample_list,fwd_master_file_list)

lower_rv_df_exp = create_and_explode_dataframe(rv_master_amplicon_number_list,rv_lower_list,
                                               rv_master_sample_list,rv_master_file_list)

lower_amplicon_df_exp = create_and_explode_dataframe(amplicon_master_number_list,amplicon_lower_list,
                                               amplicon_master_sample_list,amplicon_master_file_list)

lower_primer_amplicon_df_exp = create_and_explode_dataframe(primer_amplicon_number_list,primer_amplicon_lower_list,
                                               primer_amplicon_sample_list,primer_amplicon_master_file_list)


lower_fwd_df_exp.to_csv(flower,sep=",", index=False) # save output of weighted lowers, amplicons etc to csv
lower_rv_df_exp.to_csv(rlower,sep=",", index=False)
lower_amplicon_df_exp.to_csv(alower,sep=",", index=False)
lower_primer_amplicon_df_exp.to_csv(allcsv,sep=",",index=False)


#### plot primer_amplicon - across entire amplicon depth

plot_and_save_groups(lower_primer_amplicon_df_exp, per_plot, axhline, outdir_multi, "Primers_and_amplicon",xticks,log_condition)
plot_and_save_groups(lower_fwd_df_exp, per_plot, axhline, outdir_multi, "Forward",xticks,log_condition)
plot_and_save_groups(lower_rv_df_exp, per_plot, axhline, outdir_multi, "Reverse",xticks,log_condition)
plot_and_save_groups(lower_amplicon_df_exp, per_plot, axhline, outdir_multi, "Amplicon_only",xticks,log_condition)

print("Done :)")
