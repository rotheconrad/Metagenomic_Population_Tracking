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

import argparse
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
from scipy.stats.kde import gaussian_kde
from scipy.signal import find_peaks
import matplotlib
from matplotlib.patches import Patch
matplotlib.use('agg')
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


def func1(
    fasta, blasts, outdir, kde_bandwidth, valley_modifier,
    draw_threshold, breadth_cutoff
    ):
    "reads files builds plots"

    # Get contig lengths from reference fasta
    lengths = {}
    genome_length = 1 # total ref fasta length
    data = {
            'Sample': [], 'Location': [], 'Depth': [],
            'pident': [], 'position': [], 'rlen': [], 'alen': []
            }

    # Read through ref genome fasta and get genome length
    with open(fasta, 'r') as file:
        for name, seq in read_fasta(file):
            lengths[name[1:]] = genome_length + 1 # genome length by contig name
            genome_length += len(seq) # track total length

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
                alen = int(X[3]) # alignment length
                start = int(X[8]) # read alignment position on contig
                stop = int(X[9]) # read alignment position on contig
                pos = min(start, stop) # start position on contig
                end = max(start, stop) # stop position on contig
                rlen = int(X[12]) # length of query read
                

                # adjust alignment position by contig length
                # this is for breadth calculation and plot coordinates
                len_adjust = lengths[sid]
                adjusted_pos = pos + len_adjust - 2

                data['Sample'].append(Sample)
                data['Depth'].append(depth)
                data['Location'].append(loca)
                data['pident'].append(pid)
                data['position'].append(adjusted_pos)
                data['rlen'].append(rlen)
                data['alen'].append(alen)

    # Save data as tsv
    #mag = fasta.split('/')[-1].split('.')[0]
    df = pd.DataFrame(data)

    print('\n\nFinished processing files.\n')
    # START of Plots: write data summary tsv calculate ANIr, breadth, depth

    # start output data summary tsv
    basename = fasta.split('/')[-1].split('.')[0]
    ANIrlist = open(f'{outdir}/{basename}_stats.tsv', 'w')
    # write header for output file
    ANIrlist.write('Sample\tDetected\tDepth\tBreadth\tANIr\tMNIr\tPopCutoff\n')

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
        pop_cutoff = find_pop_threshold(subdf, samp, kde_bandwidth, valley_modifier)

        # subset population reads above population cutoff
        popx = subdf[subdf.pident >= pop_cutoff]
        if len(popx) > 0:
            # calculate depth
            depth = popx.rlen.sum() / genome_length
            # calculate breadth
            genome_array = [0] * genome_length
            # populate genome array for breadth calculation
            for alen, posi in zip(popx['alen'], popx['position']):
                for i in range(alen):
                    '''
                    print(
                        'Genome Length:', len(genome_array), 'Position:',
                        posi+i, 'Remaining Length:', alen-i
                        )
                    '''
                    try:
                        genome_array[posi+i] += 1
                    except:
                        print(f'Caution Line 160! Length Remaining: {alen-i}')
            # Calculate breadth
            breadth = sum(i > 0 for i in genome_array) / len(genome_array) * 100
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
                f'{samp}\t1\t{depth}\t{breadth}\t{anir}\t{mnir}\t{pop_cutoff}\n'
                )
            statsline = (
                f'M: {mnir:.2f}% | A: {anir:.2f}%\n'
                f'D: {depth:.2f}X | B: {breadth:.2f}%'
                )
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
        plt.savefig(f'{outdir}/{samp}_RecPlot_Mini.png', dpi=300)
        plt.close() 

    # close output file
    ANIrlist.close()


def find_pop_threshold(subdf, samp, kde_bandwidth, valley_modifier):
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
        help='Optional: Set breadth threshold detection limit (default=50.0)',
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
