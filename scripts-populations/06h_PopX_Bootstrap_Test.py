#!/usr/bin/env python

'''Sample population ANIr bootstrapping tests

Builds ANIr distribution by bootstrapping subsamples from the sequence
id of read alignments and tests probability of seeing the target
population ANIr value given the test population read alignments.

Input:

    *popx_readmap_values.tsv from 06f_TabBlast_RecPlot_Mini_Auto_v3.py

Output:

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: August 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import norm
from string import ascii_uppercase as letters
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns


def parse_input(infile): 
    
    data = {}
    pidents = {}

    with open(infile, 'r') as file:
        header = file.readline()

        for line in file:
            X = line.rstrip().split('\t')
            sample = X[0]
            Pop_Cutoff = float(X[1])
            EMR = float(X[2])
            BreadthTotal = float(X[3])
            DepthTotal = float(X[4])
            Breadth100 = float(X[5])
            Depth100 = float(X[6])
            Breadth99 = float(X[7])
            Depth99 = float(X[8])
            pids = [float(i) for i in X[9].split(',')]
            sample_size = len(pids)
            ANIr = np.mean(pids)
            MNIr = np.median(pids)

            data[sample] = {
                        'Pop_Cutoff': Pop_Cutoff,
                        'ANIr': ANIr,
                        'MNIr': MNIr,
                        'EMR': EMR,
                        'BreadthTotal': BreadthTotal,
                        'DepthTotal': DepthTotal,
                        'Breadth100': Breadth100,
                        'Depth100': Depth100,
                        'Breadth99': Breadth99,
                        'Depth99': Depth99,
                        'sample_size': sample_size,
                        }

            pidents[sample] = pids

    return data, pidents


def bootstrapping(data, pidents, bootsize, bootnumb):
    
    # initialize
    np.random.seed(42)
    bootstraps = defaultdict(list)    

    # Begin compute bootstraps for all samples in loop
    for sample, entry in data.items():

        print('\n', sample)

        # collect population measurements
        X = pidents[sample]
        n = int(entry['sample_size'] * bootsize)

        # run the bootstraps
        for i in range(bootnumb):
            if (i+1) % (bootnumb / 10) == 0:
                print(f'\t\tBootstrapping boot {i+1}')

            boot = np.random.choice(X, size=n)
            bootstraps[sample].append(np.mean(boot))

        # Calculate statistics and add to data dict
        a = np.array(bootstraps[sample], dtype=np.float64)
        n = len(a)
        x = a.sum() / n
        ss = np.abs(a - x)**2
        s = np.sqrt( ss.sum() / (n-1) )

        sim_mean = x
        sim_stdev = s
        sim_std3 = s * 3
        sim_min = np.min(a)
        sim_max = np.max(a)

        data[sample]['sim_mean'] = sim_mean
        data[sample]['sim_stdev'] = sim_stdev
        data[sample]['sim_std3'] = sim_std3
        data[sample]['sim_min'] = sim_min
        data[sample]['sim_max'] = sim_max

    return data, bootstraps


def build_plots(bootstraps, data, target_sample, outpre):

    # Pull out the target sample data and set it up
    target = data.pop(target_sample)
    target_mean = target['sim_mean']
    target_stdev = target['sim_stdev']
    target_std3 = target['sim_std3']
    target_min = target['sim_min']
    target_max = target['sim_max']

    # Begin Plotting loop of all samples
    for sample, entry in data.items():

        print('\n', 'Plotting sample:', sample)

        # define params for the sample
        sample_mean = entry['sim_mean']
        sample_stdev = entry['sim_stdev']
        sample_std3 = entry['sim_std3']
        sample_min = entry['sim_min']
        sample_max = entry['sim_max']

        # Compute normal distribution and PDF for sample
        sample_norm = norm(sample_mean, sample_stdev)
        sample_range = [i for i in np.arange(sample_min, sample_max, 0.001)]
        sample_pdf = [sample_norm.pdf(value) for value in sample_range]

        # Calculate probability of seeing sample mean given target mean
        if target_mean >= sample_mean:
            std3 = sample_mean + sample_std3
            pvalue = 1 - sample_norm.cdf(target_mean)
            htextalign = 'right'
        else:
            std3 = sample_mean - sample_std3
            pvalue = sample_norm.cdf(target_mean)
            htextalign = 'left'

        # add pvalue to data dict
        data[sample]['sim_pvalue'] = pvalue

        # define colors
        #acolor = '#252525'
        bcolor = '#969696'
        pcolor = '#000000'
        mcolor = '#525252'
        gridM = '#bdbdbd'
        alpha = 0.3

        # setup to plot
        fig, ax = plt.subplots(figsize=(8.5, 4), dpi=300)

        # plot histogram of sample bootstraps
        ax.hist(
            bootstraps[sample], bins=20, density=True, color=bcolor, alpha=alpha
            )
        # plot PDF based on sample bootstraps
        ax.plot(sample_range, sample_pdf, color=pcolor)

        # set plot/grid style
        ax.minorticks_on()
        ax.tick_params(
            which='minor', axis='both', left=False, bottom=False
            )
        ax.tick_params(
                    which='major', axis='both',
                    left=True, bottom=True,
                    size=6, width=2, tickdir='inout',
                    labelsize=11, zorder=10
                    )
        ax.yaxis.grid(
                which="major", color=gridM, linestyle='--',
                linewidth=1.5, alpha=0.2, zorder=1
                )
        ax.set_axisbelow(True)
        for spine in ax.spines.values(): spine.set_linewidth(2)

        # plot std3 and measured diff lines with pvalues
        axes = plt.gca()
        y_min, y_max = axes.get_ylim()
        ytextpos = y_max - (y_max / 20)

        ax.axvline(
            x=std3, ymin=0, ymax=1, color=mcolor, linewidth=2,
            linestyle='--', label='3 Standard Deviations'
            )
        std3_text = (
            f' ANIr at 3 Stdevs: {std3:.2f} \n'
            f' ANIr of Target: {target_mean:.2f} \n'
            f' pvalue: {pvalue:.6f} '
            )
        ax.text(
            std3, ytextpos, std3_text, color=mcolor, fontsize=11,
            horizontalalignment=htextalign, verticalalignment='top'
            )

        # set plot title and axis labels
        ax.set_title(
            f'{sample}: Bootstrapped ANIr',
            fontsize=18, fontweight='bold', y=1.02
            )
        ax.set_xlabel(
            'Bootstrapped ANIr Values',
            fontsize=14, fontweight='bold', y=-0.02
            )
        ax.set_ylabel(
            'PDF',
            fontsize=14, fontweight='bold', x=-0.02
            )

        # save figure and close
        plt.tight_layout()
        plt.savefig(f'{outpre}_{sample}_PopX_Bootstrap.png')
        plt.close() 

    return data


def write_data_file(data, target_sample, outpre):

    # setup to build and export dataframe as tsv
    df = defaultdict(list)

    # add test population data
    for sample, entry in data.items():
        df['Sample'].append(sample)
        df['Boot_Mean'].append(entry['sim_mean'])
        df['Boot_Stdev'].append(entry['sim_stdev'])
        df['Boot_Min'].append(entry['sim_min'])
        df['Boot_Max'].append(entry['sim_max'])
        df['Boot_pvalue'].append(entry['sim_pvalue'])

    # add target population placeholder
    df['Sample'].append(target_sample)
    df['Boot_Mean'].append('Target')
    df['Boot_Stdev'].append('Target')
    df['Boot_Min'].append('Target')
    df['Boot_Max'].append('Target')
    df['Boot_pvalue'].append('Target')

    df = pd.DataFrame(df)

    df.to_csv(f'{outpre}_bootstrap_data.tsv', index=False, sep='\t')

    return df


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_popx_file',
        help='Please specify the popx input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-t', '--target_sample',
        help='Please specify the sample name with the target population!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the filename prefix for the plots!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-b', '--bootstrap_boot_size',
        help='OPTIONAL: Size of sample to bootstrap \
                (Default=0.02 * target population sample size).',
        metavar='',
        type=float,
        default=0.02,
        required=False
        )
    parser.add_argument(
        '-n', '--bootstrap_number_boots',
        help='OPTIONAL: Number of boots to strap (Default=10000).',
        metavar='',
        type=int,
        default=10000,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script ...\n')

    # define input parameters
    infile = args['input_popx_file']
    target_sample = args['target_sample']
    outpre = args['output_prefix']
    bootsize = args['bootstrap_boot_size']
    bootnumb = args['bootstrap_number_boots']

    # Parse the input file
    data, pidents = parse_input(infile)
    print('\n\nFinished parsing input file. Bootstrapping ...\n')

    # bootstrap the target population
    data, bootstraps = bootstrapping(
                    data, pidents, bootsize, bootnumb
                    )
    print('\n\nFinished Bootstrapping. Plotting ...\n')

    # Build and write out the plot
    data = build_plots(bootstraps, data, target_sample, outpre)

    print('\n\nFinished plotting. Writing output tsv ...')

    # Build results dataframe and write out tsv
    df = write_data_file(data, target_sample, outpre)

    print('\n\nComplete success space cadet!! We finished without errors.\n\n')


if __name__ == "__main__":
    main()
