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
fmean = "Fmean.csv"
rmean = "Rmean.csv"
amean = "Amean.csv"
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

Find gaps in primer scheme data - returns mean depth



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
argParser.add_argument("-a", "--all_output", type=str, default='Primer_and_amplicon_depth.csv', help="csv file name for mean output for depth across primer and amplicon - default is 'Primer_and_amplicon_depth.csv'")
argParser.add_argument("-f", "--fmean_output", type=str, default='Forward_primer_depth.csv', help="csv file name for fwd primer mean output - default is 'Forward_primer_depth.csv'")
argParser.add_argument("-r", "--rmean_output", type=str, default='Reverse_primer_depth.csv', help="csv file name for rv primer mean output - default is 'Reverse_primer_depth.csv'")
argParser.add_argument("-am", "--amean_output", type=str, default='Amplicon_depth.csv', help="csv file name for amplicon mean output - default is 'Amplicon_depth.csv'")
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
fmean = args.fmean_output
rmean = args.rmean_output
amean = args.amean_output
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

def calculate_mean(sub_dfs, column_name): ## Will calculate the mean depth across each amplicon
    mean_list = []
    for sub_df in sub_dfs:
        mean = round(sub_df[column_name].mean())
        mean_list.append(mean)
    return mean_list

def create_new_df(mean_dict, primer_number): ## creates a new df with mean for each primer
    new_df = pd.DataFrame(mean_dict)
    new_df["primer_number"] = primer_number
    new_df

    return new_df

def create_and_explode_dataframe(amplicon_number_list, mean_list, sample_list, file_list): ## makes df with mean depths and formats
    df = pd.DataFrame()
    df['Amplicon_number'] = amplicon_number_list
    df['MeanDepth'] = mean_list
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
            ax.scatter(group["Amplicon_number"], group["MeanDepth"], label=name)

            if log == "log":
                plt.yscale("log", base=10)

            plt.xlabel("Amplicon number")
            plt.ylabel("Mean Depth")
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

fwd_mean_list = []
rv_mean_list = []
amplicon_mean_list = []
primer_amplicon_mean_list = []

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

    fwd_mean = calculate_mean(fwd_sub_dfs, column_name) ## calculate mean for each
    rv_mean = calculate_mean(rv_sub_dfs, column_name)
    amplicon_mean = calculate_mean(amp_sub_dfs, column_name)
    primer_amplicon_mean = calculate_mean(primer_amp_sub_dfs, column_name)

    fwd_mean_list.append(fwd_mean) ## append into list
    rv_mean_list.append(rv_mean)
    amplicon_mean_list.append(amplicon_mean)
    primer_amplicon_mean_list.append(primer_amplicon_mean)

    fwd_master_amplicon_number_list.append(amplicon_num_list)
    rv_master_amplicon_number_list.append(amplicon_num_list)
    amplicon_master_number_list.append(amplicon_num_list)
    primer_amplicon_number_list.append(amplicon_num_list)


## get exploded dfs

mean_fwd_df_exp = create_and_explode_dataframe(fwd_master_amplicon_number_list,fwd_mean_list,
                                               fwd_master_sample_list,fwd_master_file_list)

mean_rv_df_exp = create_and_explode_dataframe(rv_master_amplicon_number_list,rv_mean_list,
                                               rv_master_sample_list,rv_master_file_list)

mean_amplicon_df_exp = create_and_explode_dataframe(amplicon_master_number_list,amplicon_mean_list,
                                               amplicon_master_sample_list,amplicon_master_file_list)

mean_primer_amplicon_df_exp = create_and_explode_dataframe(primer_amplicon_number_list,primer_amplicon_mean_list,
                                               primer_amplicon_sample_list,primer_amplicon_master_file_list)


mean_fwd_df_exp.to_csv(fmean,sep=",", index=False) # save output of weighted means, amplicons etc to csv
mean_rv_df_exp.to_csv(rmean,sep=",", index=False)
mean_amplicon_df_exp.to_csv(amean,sep=",", index=False)
mean_primer_amplicon_df_exp.to_csv(allcsv,sep=",",index=False)


#### plot primer_amplicon - across entire amplicon depth

plot_and_save_groups(mean_primer_amplicon_df_exp, per_plot, axhline, outdir_multi, "Primers_and_amplicon",xticks,log_condition)
plot_and_save_groups(mean_fwd_df_exp, per_plot, axhline, outdir_multi, "Forward",xticks,log_condition)
plot_and_save_groups(mean_rv_df_exp, per_plot, axhline, outdir_multi, "Reverse",xticks,log_condition)
plot_and_save_groups(mean_amplicon_df_exp, per_plot, axhline, outdir_multi, "Amplicon_only",xticks,log_condition)

print("Done :)")
