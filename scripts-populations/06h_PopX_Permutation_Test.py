#!/usr/bin/env python

'''Sample population mean ANIr difference permutation tests

Tests the null hypothesis that the target sample population and the
test sample population are the same.

Concatenates population measurements (sequence id of read alignment) from
both samples to simulate a single population and runs permutations tests
calculating the difference of means between two samples given that they
come from the same population.

Calculates pvalue of measured ANIr difference against the simulated
ANIr difference distribution.

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


def readmap_permutation_test(data, pidents, target_sample, numbperms):

    # initialize
    np.random.seed(42)
    simulated_ANIr_diffs = defaultdict(list)

    # define target population
    target_sample_size = data[target_sample]['sample_size']
    target_pop_values = pidents[target_sample]
    target_pop_ANIr = data[target_sample]['ANIr']
    data[target_sample]['measured_diff'] = 'Target'
    data[target_sample]['perm_mean'] = 'Target'
    data[target_sample]['perm_std'] = 'Target'
    data[target_sample]['pvalue'] = 'Target'
    data[target_sample]['std3'] = 'Target'
    data[target_sample]['std3_pvalue'] = 'Target'

    # Begin compute permutations for all samples loop
    for sample, pop in data.items():
        if sample == target_sample: continue
        print('\n', sample)
        # define test population
        test_sample_size = data[sample]['sample_size']
        test_pop_values = pidents[sample]
        test_pop_ANIr = data[sample]['ANIr']

        measured_ANIr_diff = np.abs(target_pop_ANIr - test_pop_ANIr)

        combined_values = target_pop_values + test_pop_values
        z = target_sample_size

        # do the permutations
        for i in range(numbperms):
            if (i + 1) % (numbperms / 10) == 0:
                print(f'\t\tPermuations completed {i+1} ...')
            prm = list(np.random.permutation(combined_values))
            A = prm[:z]
            B = prm[z:]

            simdiff = np.mean(A) - np.mean(B)
            simulated_ANIr_diffs[sample].append(simdiff)

        # compute statistics and pvalue
        a = np.array(simulated_ANIr_diffs[sample], dtype=np.float64)
        n = len(a)
        x = a.sum() / n
        ss = np.abs(a - x)**2
        s = np.sqrt( ss.sum() / (n - 1) )

        mean_simulated_diff = x
        stdev_simulated_diff = s

        std3 = x + (s * 3)

        pvalue_measured_diff = 1 - norm.cdf(
                measured_ANIr_diff, mean_simulated_diff, stdev_simulated_diff
                )

        std3_pvalue = 1 - norm.cdf(
                std3, mean_simulated_diff, stdev_simulated_diff
                )

        # store values
        data[sample]['measured_diff'] = measured_ANIr_diff
        data[sample]['perm_mean'] = mean_simulated_diff
        data[sample]['perm_std'] = stdev_simulated_diff
        data[sample]['pvalue'] = pvalue_measured_diff
        data[sample]['std3'] = std3
        data[sample]['std3_pvalue'] = std3_pvalue

    return simulated_ANIr_diffs, data


def write_data_file(perms, data, target_sample, outpre):

    # compute median values for samples in data and setup dataframe.
    df = defaultdict(list)
    for sample, entry in data.items():
        df['Sample'].append(sample)
        df['ANIr'].append(entry['ANIr'])
        df['MNIr'].append(entry['MNIr'])
        df['EMR'].append(entry['EMR'])
        df['BreadthTotal'].append(entry['BreadthTotal'])
        df['DepthTotal'].append(entry['DepthTotal'])
        df['Breadth100'].append(entry['Breadth100'])
        df['Depth100'].append(entry['Depth100'])
        df['Breadth99'].append(entry['Breadth99'])
        df['Depth99'].append(entry['Depth99'])
        df['sample_size'].append(entry['sample_size'])
        df['Measured_Diff'].append(entry['measured_diff'])
        df['Perm_Mean_Diff'].append(entry['perm_mean'])
        df['Perm_diff_Stdev'].append(entry['perm_std'])
        df['pvalue_Measured_Diff'].append(entry['pvalue'])

    df = pd.DataFrame(df)

    df.to_csv(f'{outpre}_permutation_data.tsv', index=False, sep='\t')

    return df


def build_plots(data, perms, target_sample, outpre):

    # define colors
    #acolor = '#252525'
    bcolor = '#969696'
    lcolor = '#000000'
    mcolor = '#525252'
    gridM = '#bdbdbd'
    alpha = 0.3

    # define target population
    #target_pop_values = pidents[target_sample]
    #target_pop_ANIr = data[target_sample]['ANIr']

    for sample, entry in data.items():
        if sample == target_sample: continue

        # retrieve data for each sample
        testpop_simdiffs = perms[sample]
        testpop_mean_diff = entry['perm_mean']
        testpop_diff_std = entry['perm_std']
        testpop_measured_diff = entry['measured_diff']
        testpop_pvalue = entry['pvalue']
        testpop_std3 = entry['std3']
        testpop_std3_pvalue = entry['std3_pvalue']

        # Model normal dist for testpop permution results
        dist = norm(testpop_mean_diff, testpop_diff_std)

        # compute PDF from minimum and maximum median values in dataset
        rangemin = min(testpop_simdiffs)
        rangemax = max(testpop_simdiffs)
        value_range = [i for i in np.arange(rangemin, rangemax, 0.001)]
        probabilities = [dist.pdf(value) for value in value_range]

        # setup to plot
        fig, ax = plt.subplots(figsize=(8.5, 4), dpi=300)

        # plot histogram permution simulated differences
        ax.hist(testpop_simdiffs, bins=20, density=True, color=bcolor, alpha=alpha)
        # plot PDF based on permuations
        ax.plot(value_range, probabilities, color=lcolor)

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
        measured_diff = entry['measured_diff']
        pvalue = entry['pvalue']
        std3 = entry['std3']
        std3_pvalue = entry['std3_pvalue']

        axes = plt.gca()
        y_min, y_max = axes.get_ylim()
        ytextpos = y_max - (y_max / 20)

        ax.axvline(
            x=testpop_std3, ymin=0, ymax=1, color=mcolor, linewidth=2,
            linestyle='--', label='3 Standard Deviations'
            )
        std3_text = (
            f'Mean Diff at 3 Stdevs: {std3:.6f} \n'
            f'Measured Diff: {measured_diff:.6f} \n'
            f'pvalue: {pvalue:.6f} '
            )
        ax.text(
            testpop_std3, ytextpos, std3_text, color=mcolor, fontsize=11,
            horizontalalignment='right', verticalalignment='top'
            )

        # set plot title and axis labels
        ax.set_title(
            f'{sample}: Permutation Test',
            fontsize=18, fontweight='bold', y=1.02
            )
        ax.set_xlabel(
            f'{target_sample} - {sample} ANIr difference',
            fontsize=14, fontweight='bold', y=-0.02
            )
        ax.set_ylabel(
            'PDF',
            fontsize=14, fontweight='bold', x=-0.02
            )

        # save figure and close
        plt.tight_layout()
        plt.savefig(f'{outpre}_{sample}_permutation_test.png')
        plt.close() 


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
        '-p', '--number_of_permutations',
        help='OPTIONAL: Number of permuations to run (Default=10000).',
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
    numbperms = args['number_of_permutations']

    # Parse the input file
    data, pidents = parse_input(infile)
    print('\n\nFinished parsing input file. Running permuation tests ...\n')

    # Run permutation tests
    perms, data = readmap_permutation_test(
                                data, pidents, target_sample, numbperms
                                )
    print('\n\nFinished Permutations. Writing output tsv ...\n')

    # Build results dataframe and write out tsv
    df = write_data_file(perms, data, target_sample, outpre)

    print('\n\nFinished writing output tsv. Plotting ...')

    # Build and write out the plots.
    build_plots(data, perms, target_sample, outpre)

    print('\n\nComplete success space cadet!! We finished without errors.\n\n')


if __name__ == "__main__":
    main()
