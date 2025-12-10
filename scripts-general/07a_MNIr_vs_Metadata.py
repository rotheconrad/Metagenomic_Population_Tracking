#!/usr/bin/env python

'''Builds scatter correlation plot between MNIr values and metadata

Written for GoM Project:
In puts:
    1) Metagenome_Samples_MetaData_Density.tsv
    2) Updated RecPlot Mini Auto *stats.tsv from STEP 07 number 2.

This script will write a metadata file with only the samples having a
population passing the population cutoffs that can be used to make
a TS plot for samples tailored to specific populations.

Need to install python seawater or gsw package to calculate density:
https://pypi.org/project/seawater/
conda install -c conda-forge seawater or pip install seawater

https://pypi.org/project/gsw/
conda install -c conda-forge gsw or pip install gsw

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 02 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr as corr
import pandas as pd

def gather_metadata(metadata):
    data = pd.read_csv(metadata, sep='\t')
    data['Sample'] = data['Sample'].apply(lambda x: '_'.join(x.split('_')[:3]))
    data.set_index('Sample', inplace=True)
    data = data[['Ocean', 'Density', 'Depth', 'Temp', 'Salinity']]
    data.sort_values(by='Sample', inplace=True)
    return data

def gather_popstats(popstats, metadata):
    data = pd.read_csv(popstats, sep='\t').set_index('Sample')
    data.sort_values(by='Sample', inplace=True)
    data.drop(['Depth', 'Breadth'], 1, inplace=True)
    data = data.join(metadata, how='outer').drop('PopCutoff', 1)
    data = data[data['Detected'] == 1]
    data.drop('Detected', 1, inplace=True)
    return data

def gather_correlation(data):
    # Compute Pearson Correlation Coefficient
    pdepth = corr(data.Depth, data.MNIr)
    ptemp = corr(data.Temp, data.MNIr)
    psal = corr(data.Salinity, data.MNIr)
    pdens = corr(data.Density, data.MNIr)
    df_corr = {
            'Depth': pdepth, 'Temp': ptemp, 'Salinity': psal, 'Density': pdens
            }
    return df_corr

def generate_plots(data, corr, outpre, metacolors):

    mnir = data['MNIr']
    labels = {
            'Depth': 'Depth (meters)',
            'Temp': 'Temperature (celcius)',
            'Salinity': 'Salinity (ppt)',
            'Density': 'Density (kg/m3)'
            }

    if metacolors:
        # Load the metadata
        metac = pd.read_csv(metacolors, sep='\t', header=0, dtype = str)
        colors = dict(zip(metac['Labels'], metac['Colors']))
        mappedcolor = data['Ocean'].map(colors)

    for param in ['Depth', 'Temp', 'Density', 'Salinity']:
        pearsonr, pvalue = corr[param]

        stats_line = (
            f"Pearson r: {pearsonr:.4f}\n"
            f"p-value: {pvalue:.4f}"
            )

        # Set Colors
        full_color = "#a50f15"
        full_line_color = "#000000"
        rep_color = "#000000"

        # Build the plot
        fig, ax = plt.subplots(figsize=(12, 10))

        # plot title, labels, and text
        ax.set_title(
            f'MNIr vs {labels[param]}',
            fontsize=50, y=1.02, color=rep_color
            )
        bbox_props = dict(boxstyle="round", alpha=0.5, fc='w', ec='0.5')
        ax.text(
            0.28, 0.98, stats_line,
            fontsize=24, color='#252525',
            verticalalignment='top', horizontalalignment='right',
            transform=ax.transAxes, bbox=bbox_props, size=20
            )
        ax.set_xlabel(
            f'{labels[param]}',
            fontsize=30, fontweight='bold', y=-0.02
            )
        ax.set_ylabel(
            'Median Nucleotide Identity\nof mapped reads (MNIr)',
            fontsize=30, fontweight='bold', x=-0.02
            )

        if metacolors:
            # plot the data colored by Ocean Category
            ax.scatter(
                data[param], mnir, c=mappedcolor,
                marker='o', s=300, alpha=0.50
                )

        else:
            # plot the data for df_full
            ax.plot(
                data[param], mnir, linestyle='None',
                marker='o', markersize=15, color=full_color, alpha=0.50
                )

        # set the axis parameters / style
        ax.minorticks_on()
        ax.tick_params(axis='both', labelsize=18)
        ax.tick_params(
            axis='x', which='major', direction='in', color='k',
            width=6, length=12, bottom=True, zorder=3
            )

        # set grid style
        ax.yaxis.grid(which="minor", color='#d9d9d9', linestyle='--', linewidth=1)
        ax.xaxis.grid(which="minor", color='#f0f0f0', linestyle='-', linewidth=1)
        ax.yaxis.grid(which="major", color='#d9d9d9', linestyle='--', linewidth=1.5)
        ax.xaxis.grid(which="major", color='#f0f0f0', linestyle='-', linewidth=2)
        ax.set_axisbelow(True)

        # adjust layout, save, and close
        fig.set_tight_layout(True)
        plt.savefig(f'{outpre}_{param}_corr_plot.png')
        plt.close()

    # Write separate legend
    fig, ax = plt.subplots(figsize=(10,10))
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    for label, color in colors.items():
        ax.bar(0, 0, color=color, label=label, linewidth=0, alpha=0.5)

    ax.legend(
        title='Ocean', title_fontsize='xx-large', loc="center", frameon=False,
        markerscale=5, fontsize='xx-large'
        )

    plt.savefig(f'{outpre}_Legend.png', dpi=300)
    plt.close()

    # Write metadata file with samples where population is detected.
    data.to_csv(f'{outpre}_Curated_Metadata.tsv', sep='\t')

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-m', '--metadata_input_file',
        help='Please specify the metadata input file name!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-a', '--pop_stats_file',
        help='Please specify the population stats file name!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify the output file name prefix!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--metacolors',
        help='OPTIONAL: Metadata to color samples by!',
        metavar='',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # define input params
    metadata = args['metadata_input_file']
    popstats = args['pop_stats_file']
    outpre = args['output_file_prefix']
    metacolors = args['metacolors']

    # Collect metadata
    mdf = gather_metadata(metadata)
    # Pair pop stats with metadata
    adf = gather_popstats(popstats, mdf)
    # compute correlations
    cdf = gather_correlation(adf)
    # generate plots
    generate_plots(adf, cdf, outpre, metacolors)


if __name__ == "__main__":
    main()
