#!/usr/bin/env python

'''Plots from folder of filtered tabular blasts.

Builds boxplots and hexplots.

Naming scheme for PAC and GoM metagenome samples and GoM MAGs.

For different projects change lines:
Sample, depth, and loca lines 96, 97, 98

## To merge all png files into a single pdf file in bash ##
convert directory/*png Filename.pdf

For Mac OS X: Install HomeBrew, Install ImageMagick (includes convert)
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install imagemagick

see: 
https://apple.stackexchange.com/questions/335635/where-is-the-convert-command-in-macos

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

import argparse, copy
from os import listdir
from os.path import isfile, join
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats.kde import gaussian_kde
from scipy.signal import find_peaks
import matplotlib
from matplotlib.patches import Patch
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns

# Getting false positive warning on line 134
pd.options.mode.chained_assignment = None  # default='warn'

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def truncate(x, tad=0.8):
    """ returns tad range of a list/array """

    xsorted = sorted(x)
    xlen = len(x) # get length of the list
    inverse_tad = round((1.0 - tad) / 2.0, 2) # to get top and bottom 
    q = int(xlen * inverse_tad) # to get top and bottom
    bottom = q
    top = xlen - q
    tad_range = xsorted[bottom:top] # slice list

    return tad_range


def get_average(l, tad=0.8):
    """ returns average from truncated list """

    trunc_val = truncate(l, tad)

    if sum(trunc_val) == 0:
        average = 0
    else:
        average = sum(trunc_val) / len(trunc_val)

    return average


def calc_coverage(popx, contig_dicts):
    working_dict = copy.deepcopy(contig_dicts)
    # populate genome array for breadth calculation
    for index, row in popx.iterrows():
        for i in range(row['origpos'], row['origpos']+row['poslen']+1, 1):
            working_dict[row['contig_name']][i] += 1

    # Calculate breadth
    genomecov = []
    for contig, positions in working_dict.items():
        values = list(positions.values())
        genomecov.extend(values)

    breadth = sum(i > 0 for i in genomecov) / len(genomecov) * 100
    depth = get_average(genomecov)

    return breadth, depth


def func1(
    fasta, blasts, outdir, kde_bandwidth, valley_modifier,
    draw_threshold, breadth_cutoff
    ):
    "reads files builds plots"

    # Get contig lengths from reference fasta
    contig_dicts = defaultdict(dict)
    contig_lengths = {}
    adjusted_lengths = {}
    genome_length = 0 # total ref fasta length
    data = defaultdict(list)

    # Read through ref genome fasta and get genome length
    with open(fasta, 'r') as file:
        for name, seq in read_fasta(file):
            contig_name = name.split(' ')[0][1:]
            length = len(seq) # calculate length of contig

            contig_lengths[contig_name] = length
            adjusted_lengths[contig_name] = genome_length
            genome_length += length

            # This populates the dictionary with value of zero for each
            # base pair position for each contig in the genome fasta
            for i in range(1, length+1, 1):
                contig_dicts[contig_name][i] = 0

    # Grab list of blast files from blast directory
    file_list = [f for f in listdir(blasts) if isfile(join(blasts, f))]
    # remove .DS_Store file for stupid MAC OS
    if '.DS_Store' in file_list: file_list.remove('.DS_Store')

    # Keep a sample list
    Sample_List = []
    # For Ocean depth profile samples keep track of sample depth
    # this could modified to track sample location or timepoint etc.
    depths = {}

    # Read through blast files and populate the data dict
    for file in file_list:
        print(f'Processing File: {file} ...')
        # Retreive sample data from file name
        Sample = file.split('.')[0].split('-')[0]
        depth = Sample.split('_')[2]
        loca = '_'.join(Sample.split('_')[:2])
        Sample_List.append(Sample)
        depths[Sample] = depth

        with open(f'{blasts}/{file}', 'r') as f:
            # header = f.readline() # tabular blasts shouldn't have a header!
            for line in f:
                X = line.rstrip().split('\t')
                sid = X[1] # contig name query read maps to
                pid = float(X[2]) # percent sequence identity of alignment
                #alen = int(X[3]) # alignment length
                start = int(X[8]) # read alignment position on contig
                stop = int(X[9]) # read alignment position on contig
                pos = min(start, stop) # start position on contig
                end = max(start, stop) # stop position on contig
                poslen = end - pos
                #rlen = int(X[12]) # length of query read
                
                # adjust alignment position by contig length
                # this is for breadth calculation and plot coordinates
                len_adjust = adjusted_lengths[sid]
                adjusted_pos = pos + len_adjust

                data['contig_name'].append(sid)
                data['Sample'].append(Sample)
                data['Depth'].append(depth)
                data['Location'].append(loca)
                data['pident'].append(pid)
                data['position'].append(adjusted_pos)
                data['poslen'].append(poslen)
                data['origpos'].append(pos) #unadjusted orig position
                #data['rlen'].append(rlen)
                #data['alen'].append(alen)

    # Save data as tsv
    #mag = fasta.split('/')[-1].split('.')[0]
    df = pd.DataFrame(data)

    print('\n\nFinished processing files.\n')
    # START of Plots: write data summary tsv calculate ANIr, breadth, depth

    # start output data summary tsv
    # create output directory if doesn't exist
    Path(outdir).mkdir(parents=True, exist_ok=True)
    # grab the basename of the input file
    basename = fasta.split('/')[-1].split('.')[0]
    basefile = f'{outdir}/{basename}'
    ANIrlist = open(f'{basefile}_stats.tsv', 'w')
    popx_readmap_value_file = open(f'{basefile}_popx_readmap_values.tsv', 'w')
    # write header for output file
    ANIrlist.write('Sample\tDetected\tDepth\tBreadth\tANIr\tMNIr\tPopCutoff\n')
    popx_readmap_value_file.write(
        'Sample\tPopCutoff\tExactMatchRatio(EMR)\tBreadthTotal\tDepthTotal\t'
        'Breadth100\tDepth100\tBreadth99\tDepth99\t'
        'pident above population threshold\n'
        )
    # start plot
    sns.set_style("white")
    for samp in Sample_List:
        print(f'Plotting Sample: {samp}')
        # subsample dataframe to current sample
        subdf = df[df['Sample'] == samp]
        # retrieve sample depth (or location or time etc)
        sample_depth = depths[samp]

        # modify genome positions for miniplot x-axis in Mbps
        subdf['pos'] = subdf['position'].div(1000000)
        xmax = genome_length / 1000000
        xmid = xmax / 2

        # Find the population cutoff threshold
        pop_cutoff = find_pop_threshold(
            subdf, samp, kde_bandwidth, valley_modifier, f'{outdir}/{samp}'
            )

        # subset population reads above population cutoff
        popx = subdf[subdf.pident >= pop_cutoff]
        if len(popx) > 0:
            # calculate coverage
            #depth = popx.poslen.sum() / genome_length
            breadth, depth = calc_coverage(popx, contig_dicts)

        else:
            depth = 0
            breadth = 0
            alpha = 0.0

        if breadth >= breadth_cutoff:
            # calcualte anir
            anir = popx.pident.mean()
            mnir = popx.pident.median()
            ytextpos = anir - 4
            ANIrlist.write(
                f'{samp}\t1\t{depth:.2f}\t{breadth:.2f}\t{anir:.2f}\t'
                f'{mnir:.2f}\t{pop_cutoff:.2f}\n'
                )
            statsline = (
                f'M: {mnir:.2f}% | A: {anir:.2f}%\n'
                f'D: {depth:.2f}X | B: {breadth:.2f}%'
                )
            ####################################################################
            # write data for 06h PopX Eval

            # Exact Match Ratio
            pidents = popx.pident.tolist()
            exactmatch = [i for i in pidents if i >= 99]
            EMR = round(len(exactmatch) / len(pidents), 4)

            # breadth 100
            popx100 = popx[popx.pident == 100]
            breadth100, depth100 = calc_coverage(popx100, contig_dicts)

            # breadth 99
            popx99 = popx[popx.pident >= 99]
            breadth99, depth99 = calc_coverage(popx99, contig_dicts)

            pident_string = ','.join([str(i) for i in pidents])
            popx_readmap_values = (
                f'{samp}\t{pop_cutoff:.2f}\t{EMR:.2f}\t{breadth}\t{depth}\t'
                f'{breadth100:.2f}\t{depth100:.2f}\t{breadth99:.2f}\t'
                f'{depth99:.2f}\t{pident_string}\n'
                )
            popx_readmap_value_file.write(popx_readmap_values)

        else:
            anir = 0
            mnir = 0
            ytextpos = 97
            ANIrlist.write(f'{samp}\t0\t0\t0\t0\t0\t0\n')
            statsline = '-Below Detection-'

        # define colors and limits
        hcol = '#252525'
        lcol = '#b2182b' if breadth > breadth_cutoff else '#2166ac' # '#d7191c'

        # build plot
        h = sns.JointGrid()
        _ = sns.histplot(
                    data=subdf, x="pos", y="pident",
                    bins=30, ax=h.ax_joint, cmap='Greys'
                    )
        if breadth >= breadth_cutoff:
            _ = sns.histplot(
                    data=popx, x="pos", bins=30, ax=h.ax_marg_x, color=hcol
                    )
        else:
            _ = sns.histplot(
                    data=subdf, x="pos", bins=30, ax=h.ax_marg_x, color='#FFFFFF'
                    )

        bns = [i for i in range(70,101)]
        _ = sns.histplot(
                    data=subdf, y="pident", bins=bns, ax=h.ax_marg_y, color=hcol
                    )


        h.ax_joint.axhline(y=anir, ls='--', color=lcol, lw=3)
        h.ax_joint.axhline(y=mnir, ls=':', color=lcol, lw=2)

        if draw_threshold:
            h.ax_joint.axhline(y=pop_cutoff, ls='--', color='#2166ac', lw=1)

        h.ax_joint.text(
                #xmax,
                xmid,
                ytextpos,
                statsline,
                fontsize=12, color=lcol, horizontalalignment='center',
                verticalalignment='top'
                )
        bbox_props = dict(boxstyle="round", alpha=0.7, fc='w', ec='0.5')
        h.ax_joint.text(
                    xmax-0.05, 72, sample_depth,
                    ha='right', va='bottom', bbox=bbox_props, size=14
                    )
        h.ax_joint.set_xlabel('')
        h.ax_joint.set_ylabel('')
        h.ax_joint.set_xlim([0,xmax])
        h.ax_joint.set_ylim([70,100])

        h.fig.set_figwidth(3)
        h.fig.set_figheight(2.5)
        plt.subplots_adjust(left=0.15, bottom=0.1, right=0.995, top=0.995)
        plt.savefig(f'{outdir}/{samp}_RecPlot_Mini.pdf', dpi=300)
        plt.close() 

    # close output file
    ANIrlist.close()


def find_pop_threshold(subdf, samp, kde_bandwidth, valley_modifier, basefile):
    pid = subdf['pident'].to_numpy()
    kde = gaussian_kde(pid, bw_method=kde_bandwidth)
    x = np.linspace(70, 100, 1000)
    y = kde(x)
    yn = y * -1

    #ypn = find_peaks(yn, width=50) # find the negative peak (ie the valley)
    ypn = find_peaks(yn)
    if len(ypn[0]) == 1:
        valley = x[ypn[0]][0] + valley_modifier
    elif len(ypn[0]) > 1:
        valley = x[max(ypn[0])] + valley_modifier
    elif len(ypn[0]) == 0:
        yp = find_peaks(y)
        peak = x[yp[0][0]]
        valley = peak - 5 if peak > 85 else peak + 10
    else:
        print(f'\n\nFind peaks error!!\n{samp}\n\n')
        valley = 0

    if valley > 95: valley = 95

    fig, ax = plt.subplots(figsize=(8.5, 4))
    _ = ax.hist(pid, bins=30, density=True)
    _ = ax.plot(x, y, 'r-', label='non-parametric KDE')
    _ = ax.plot(x, yn, 'k--', label='inverse non-parametric KDE')
    _ = ax.axvline(x=valley, color='b', linestyle='--')
    ax.set_title(
        f'{samp}: Population Cutoff',
        fontsize=18, fontweight='bold', y=1.02
        )
    ax.set_xlabel(
        'Percent Sequence Identity',
        fontsize=14, fontweight='bold', y=-0.02
        )
    ax.set_ylabel(
        'Kernel Density Estimate',
        fontsize=14, fontweight='bold', x=-0.02
        )
    ax.legend(loc='best')
    plt.savefig(f'{basefile}_Pop_Cutoff.png', dpi=300)
    plt.close()

    return valley


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-f', '--input_fasta_file',
        help='Please specify the reference fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-b', '--input_blast_dir',
        help='Please specify the blast directory!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_dir',
        help='Please specify the output directory!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-k', '--kde_bandwidth',
        help='Optional: Change bandwidth to control smoothing (default=0.25).',
        metavar='',
        type=float,
        required=False,
        default=0.3
        )
    parser.add_argument(
        '-v', '--valley_modifier',
        help='Optional: Nudge the cutoff line (default=3.0).',
        metavar='',
        type=float,
        required=False,
        default=3.0
        )
    parser.add_argument(
        '-d', '--draw_threshold',
        help='Use this flag to draw population cutoff threshold on the plot.',
        required=False,
        action='store_true'
        )
    parser.add_argument(
        '-t', '--breadth_threshold',
        help='Optional: Set breadth threshold detection limit (default=10.0)',
        metavar='',
        type=float,
        required=False,
        default=10.0
        )
    args=vars(parser.parse_args())
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script ...\n')
    func1(
        args['input_fasta_file'],
        args['input_blast_dir'],
        args['out_dir'],
        args['kde_bandwidth'],
        args['valley_modifier'],
        args['draw_threshold'],
        args['breadth_threshold']
        )


if __name__ == "__main__":
    main()
