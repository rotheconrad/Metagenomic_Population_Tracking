#!/usr/bin/env python

'''Build heatmap of MAG Biogeography

Takes a directory of *stats.tsv files as input generated from:
06f_TabBlast_RecPlot_Mini_Auto_v2.py

Writes out heatmaps for Presence, Depth, Breadth, ANIr, and MNIr.

### Future improvement. Take metadata file with samples by category
and a second file with the colors to use for each category.
Map colors to each category and write out the legends.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Jan 04 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, re
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.patches import Patch
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Set the default sans-serif font to Helvetica
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# Set default font to sans-serif
matplotlib.rcParams['font.family'] = "sans-serif"
# If you don't have Helvetica and need it see this:
# https://fowlerlab.org/2019/01/03/changing-the-sans-serif-font-to-helvetica/
# If you don't need Helvetica font just delete the above code.

def collect_stats(
    input, breadth_threshold, depth_threshold, sample_threshold, outpre, sorder
    ):

    # Grab list of blast files from blast directory
    file_list = [f for f in listdir(input) if isfile(join(input, f))]
    # remove .DS_Store file for stupid MAC OS
    if '.DS_Store' in file_list: file_list.remove('.DS_Store')

    # use counter to set sample names at first file
    count = 0

    # Initialize dataframe for each variable:
    detected = {}
    depth = {}
    breadth = {}
    anir = {}
    mnir = {}

    # Substitute numbers for MAG Names
    MAG_Names = {
            'MAG_Name': [], 'Pass_Thresholds': [], 'Sample_Count': [],
            #'MAG_Number': [], 
            }

    # Read through stats files and build dataframes of MAG x Sample
    # for each parameter Detected, Depth, Breadth, ANIr, and MNIr.
    for i, file in enumerate(file_list):
        # track mag number and substitute mag number for name
        mag_name = '_'.join(file.split('_')[:5])
        #mag = f'MAG-{i}'
        MAG_Names['MAG_Name'].append(mag_name)
        #MAG_Names['MAG_Number'].append(mag)

        # read in the file and sort it
        df = pd.read_csv(f'{input}/{file}', sep='\t', header=0)
        df.sort_values(by='Sample', inplace=True)

        # check thresholds
        b_thresh = sum(df['Breadth'] >= breadth_threshold) >= sample_threshold
        d_thresh = sum(df['Depth'] >= depth_threshold) >= sample_threshold
        s_thresh = df['Detected'].sum() >= sample_threshold

        thresholds = [b_thresh, d_thresh, s_thresh]

        # check if mag passes thresholds.
        if sum(thresholds) < 3:
            #print(f'{mag}  !! Did not pass thresholds !!')
            MAG_Names['Pass_Thresholds'].append(0)
            MAG_Names['Sample_Count'].append(0)

        # if pass threshholds add mag to data
        else:

            if count == 0:
                # Add sample names as first column
                detected['Sample'] = df['Sample'].tolist()
                depth['Sample'] = df['Sample'].tolist()
                breadth['Sample'] = df['Sample'].tolist()
                anir['Sample'] = df['Sample'].tolist()
                mnir['Sample'] = df['Sample'].tolist()
                count += 1

            # Add data for mag to dict
            detected[mag_name] = df['Detected'].tolist()
            depth[mag_name] = df['Depth'].tolist()
            breadth[mag_name] = df['Breadth'].tolist()
            anir[mag_name] = df['ANIr'].tolist()
            mnir[mag_name] = df['MNIr'].tolist()

            # Update records for summary
            MAG_Names['Pass_Thresholds'].append(1)
            MAG_Names['Sample_Count'].append(df['Detected'].sum())

    # Convert to dataframes and write tsv output files
    DTCTD = pd.DataFrame(detected).set_index('Sample')
    DPTH = pd.DataFrame(depth).set_index('Sample')
    BRDTH = pd.DataFrame(breadth).set_index('Sample')
    ANIR = pd.DataFrame(anir).set_index('Sample')
    MNIR = pd.DataFrame(mnir).set_index('Sample')
    
    # Reorder dataframes if sorder true
    if sorder:
        DTCTD = DTCTD.rename(index=sorder).sort_index()
        DPTH = DPTH.rename(index=sorder).sort_index()
        BRDTH = BRDTH.rename(index=sorder).sort_index()
        ANIR = ANIR.rename(index=sorder).sort_index()
        MNIR = MNIR.rename(index=sorder).sort_index()

    # Write the dataframes to tsv files
    DTCTD.to_csv(f'{outpre}_Presence_Absence.tsv', sep='\t')
    DPTH.to_csv(f'{outpre}_Sequencing_Depth.tsv', sep='\t')
    BRDTH.to_csv(f'{outpre}_Sequencing_Breadth.tsv', sep='\t')
    ANIR.to_csv(f'{outpre}_Population_ANIr.tsv', sep='\t')
    MNIR.to_csv(f'{outpre}_Population_MNIr.tsv', sep='\t')

    # Sort the MAG Key descending by 'Sample_Count' and write to file
    mag_key = pd.DataFrame(MAG_Names)
    mag_key.sort_values(by='Sample_Count', ascending=False, inplace=True)
    mag_key.to_csv(f'{outpre}_MAG_SampleCount.tsv', sep='\t', index=None)

    # Filtering Threshold Summary Report:
    total = len(MAG_Names['Pass_Thresholds'])
    passed = sum(MAG_Names['Pass_Thresholds'])
    failed = total - passed

    print(
        f'\n\nSummary of threshold filtering results:\n\n'
        f'\t\tMAGs need >= {breadth_threshold}% Sequencing Breadth\n'
        f'\t\tMAGs need >= {depth_threshold}X Sequencing Depth\n'
        f'\t\tin >= {sample_threshold} metagenome samples.\n\n'
        f'\t\tTotal Passed: {passed} / {total} ({passed/total*100:.2f}%)\n'
        f'\t\tTotal Failed: {failed} / {total} ({failed/total*100:.2f}%)'
        )

    # Return a dict of the dataframes for each parameter.
    return {'A': DTCTD, 'B': DPTH, 'C': BRDTH, 'D': ANIR, 'E': MNIR}


def build_heatmaps(data, metadata, metacolors, outpre, sorder, classifications):

    labels = {
            'A': ['Presence_Absence', ["#ffffff", "#000000"]],
            'B': ['Sequencing_Depth', 'Greys'],
            'C': ['Sequencing_Breadth', 'Greys'],
            'D': ['Population_ANIr', 'Greys'],
            'E': ['Population_MNIr', 'Greys']
    }

    # if metadata and colors read in and prepare the data.
    if metadata and metacolors:
        metad = pd.read_csv(metadata, sep='\t', header=0, dtype = str)
        metad.sort_values(by='Sample', inplace=True)
        # The indexes need to match metad to df.
        metad.set_index('Sample', inplace=True)
        if sorder: metad = metad.rename(index=sorder)
        # make a dict from Labels: Colors and map to metad values.
        metac = pd.read_csv(metacolors, sep='\t', header=0, dtype = str)
        metac_dict = dict(zip(metac['Labels'], metac['Colors']))
        mappedmeta = metad.stack().map(metac_dict).unstack()

    cbarpos = (0.05, .9, .2, .05) #(left, bottom, width, height),
    plotsize = (25,12)

    for key, df in data.items():
        label = labels[key][0]
        color = labels[key][1]

        print(f'\t\t\tPlotting: {label} ...')


        # set low and high cbar values
        if label == 'Presence_Absence':
            low = 0
            high = 1
        else:
            dflist = [i for i in df.to_numpy().flatten() if i != 0]
            low = min(dflist)
            high = max(dflist)
            print(f'\t\t\t\tlow: {low}\n\t\t\t\thigh: {high}')

        # plot clustermap
        if metadata and metacolors:
            sns.clustermap(
                df, figsize=plotsize, cmap=color, vmin=low, vmax=high,
                cbar_pos=cbarpos, row_colors=mappedmeta,
                cbar_kws={"orientation": "horizontal"}
                )
        else:
            sns.clustermap(
                df, figsize=plotsize, cmap=color,
                cbar_pos=cbarpos, vmin=low, vmax=high,
                cbar_kws={"orientation": "horizontal"}
                )

        plt.savefig(f'{outpre}_{label}_clustermap.pdf', dpi=300)
        plt.close()

        # plot heatmap
        if metadata and metacolors:
            #sns.set(font_scale=0.5)
            g= sns.clustermap(
                df, figsize=plotsize, cmap=color, vmin=low, vmax=high,
                cbar_pos=cbarpos, row_colors=mappedmeta,
                cbar_kws={"orientation": "horizontal"},
                row_cluster=False, xticklabels=True
                )
        else:
            sns.clustermap(
                df, figsize=plotsize, cmap=color,
                cbar_pos=cbarpos, vmin=low, vmax=high,
                cbar_kws={"orientation": "horizontal"},
                row_cluster=False
                )

        # Print Mag numbers
        plt.savefig(f'{outpre}_{label}_heatmap.pdf', dpi=300)
        plt.close()

    xlabo = [i.get_text() for i in g.ax_heatmap.get_xticklabels()]
    switch = [classifications[i] for i in xlabo]
    with open(f'{outpre}_Classifications.txt', 'w') as clssfy_out:
            clssfy_out.writelines(f'{clss}\n' for clss in switch)

    ## plot meta data legends #########################################
    if metadata and metacolors:
        for meta in metad.columns:
            labels = natural_sort(metad[meta].unique())

            fig, ax = plt.subplots(figsize=(10,10))
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

            for label in labels:
                ax.bar(
                    0, 0, color=metac_dict[label], label=label, linewidth=0
                    )

            ax.legend(
                title=meta, title_fontsize='xx-large', loc="center",
                frameon=False, markerscale=5, fontsize='xx-large'
                )

            plt.savefig(f'{outpre}_Legend_{meta}.pdf', dpi=300)
            plt.close()
    ###################################################################

def natural_sort(l): 
    # copied from:
    # https://stackoverflow.com/questions/11150239/python-natural-sorting
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
    

def get_order(samp_order):
    # read in sample order list
    sorder = {}

    with open(samp_order, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            samp = X[0]
            order = int(X[-1])
            sorder[samp] = order

    return sorder

def get_classifications(infile):
    # parse MAG Classification file
    classifications = {}

    with open(infile, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            MAG = X[0]
            classify = ' '.join(X[1].split('_')[:2])
            AAI = float(X[2])
            classifications[MAG] = f'{MAG}\t{classify}\t({AAI:.2f})'

    return classifications

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_stats_dir',
        help='Please specify the directory with *stats.tsv filse!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-md', '--metadata',
        help='OPTIONAL: Metadata to color samples by! requires -mc',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-mc', '--metacolors',
        help='OPTIONAL: Metadata to color samples by! requires -md',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the filename prefix for the plots!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-b', '--breadth_threshold',
        help='Optional: Set breadth threshold detection limit (default>=10.0)',
        metavar='',
        type=float,
        required=False,
        default=10.0
        )
    parser.add_argument(
        '-d', '--depth_threshold',
        help='Optional: Set depth threshold detection limit (default>=0.01)',
        metavar='',
        type=float,
        required=False,
        default=0.01
        )
    parser.add_argument(
        '-s', '--sample_threshold',
        help='Optional: Set sample threshold detection limit (default>=2)',
        metavar='',
        type=int,
        required=False,
        default=2
        )
    parser.add_argument(
        '-r', '--sample_reorder',
        help='Optional: Specify the order of samples in the heatmap',
        metavar='',
        type=str,
        required=False,
        )
    parser.add_argument(
        '-c', '--MAG_classification',
        help='Optional: Specify the Classification for MAGs',
        metavar='',
        type=str,
        required=False,
        )
    args=vars(parser.parse_args())
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script ...\n')

    # Check for sample order file and process if true
    sorder = None
    classifications = None
    if args['sample_reorder']: sorder = get_order(args['sample_reorder'])

    if args['MAG_classification']:
        classifications = get_classifications(args['MAG_classification'])

    # Collect stats data
    data = collect_stats(
                        args['input_stats_dir'],
                        args['breadth_threshold'],
                        args['depth_threshold'],
                        args['sample_threshold'],
                        args['output_prefix'],
                        sorder
                        )
    print('\n\nFinished processing stats files. Plotting ...\n')

    # Plot stats data
    build_heatmaps(
                data,
                args['metadata'],
                args['metacolors'],
                args['output_prefix'],
                sorder,
                classifications
                )

    print('\n\nComplete success space cadet!! We finished without errors.\n\n')


if __name__ == "__main__":
    main()
