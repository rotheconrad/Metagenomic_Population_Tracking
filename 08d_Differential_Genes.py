#!/usr/bin/env python

'''Differential gene coverage between samples.

This script works downstream of 08b_BlastPlus_CoverageMagic_Basic.py.

It reads the output directory and parses the *_gene_tad.tsv files.
It calculates the difference in normalized tad value between the
target sample population and the test sample population.

It writes out a histogram of these values with the top and bottom 2.5%
percentile marked.

It writes a tsv file containing the genes in the top and bottom 2.5%
percentile of differential gene coverage.

* All input files needs to use the same basenames.

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
from os import listdir
from os.path import isfile, join
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


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


def parse_annotations(infile):

    annotations = {}

    with open(infile, 'r') as file:
        header = file.readline()
        H = header.rstrip().split('\t')
        annotations['Header'] = '\t'.join(H[1:])

        for line in file:
            X = line.rstrip().split('\t')
            gene = X[0]
            info = '\t'.join(X[1:])
            annotations[gene] = info

    return annotations


def parse_fasta(infile):
    
    genefasta = {}

    with open(infile, 'r') as file:
        for name, seq in read_fasta(file):
            gene = name[1:].split(' ')[0]
            genefasta[gene] = seq

    return genefasta


def parse_input(indir):

    # initialize data dict:
    data = {} # lognormtads (lgnmtd)
    data2 = {} # normtads     (nmtd)
    data3 = {} # tads           (td)

    # Grab list of input files from input directory
    ext = 'gene_tad.tsv'
    file_list = [f for f in listdir(indir) if isfile(join(indir, f))]
    file_list = [f for f in file_list if '_'.join(f.split('_')[-2:]) == ext]
    # remove .DS_Store file for stupid MAC OS
    if '.DS_Store' in file_list: file_list.remove('.DS_Store')

    # Begin loop: Read through files and populate data dict
    for file in file_list:
        print(f'\t\tParsing file: {file} ...')
        lognormtads = {} # initialize dict for genes
        normtads = {}
        tads = {}

        Sample = file.split('-')[0] # define sample name
        # loop through each line in file and retrieve the normalized tad
        with open(f'{indir}/{file}', 'r') as f:
            header = f.readline() # dump the header
            for line in f:
                X = line.rstrip().split('\t')
                gene = X[0]
                tad = X[1]
                # Gene TAD normalized by genome average in coverage magic script
                # Here we take the log of that ratio to even out the difference
                # between 0 to 1 and 1 to infinity that occurs with ratios.
                norm_tad = float(X[2])
                normtads[gene] = norm_tad
                tads[gene] = tad

        # Find minimum norm tad value to adjust zeros for log ratio
        n = min([i for i in normtads.values() if i != 0])/2
        for gene, norm_tad in normtads.items():
            if norm_tad == 0: norm_tad = n
            log_norm_tad = np.log(norm_tad)
            lognormtads[gene] = log_norm_tad

        # Add genetad to data dict by sample name
        data[Sample] = lognormtads
        data2[Sample] = normtads
        data3[Sample] = tads

    return data, data2, data3 # (lgnmtd, nmtd, td)


def analyze_data(genetads, sigval, target_sample):
    
    # initialize data dict:
    data = {}
    genediffs = {}

    # grab target sample
    target = genetads.pop(target_sample)

    # read through data for each sample
    for sample, genes in genetads.items():
        stats = {}
        # Calculate difference per gene target - current
        #diffs = [target[gene] - tad for gene, tad in genes.items()]
        diff_array = []
        diff_dict = {}

        for gene, tad in genes.items():
            x = target[gene] - tad
            diff_array.append(x)
            diff_dict[gene] = x

        # Compute high and low percentiles
        #hq = round(np.quantile(diffs, (1-sigval)), 2) # high
        #lq = round(np.quantile(diffs, sigval), 2) # low
        hq = round(np.quantile(diff_array, (1-sigval)), 2) # high
        lq = round(np.quantile(diff_array, sigval), 2) # low
        # how man genes > hq and < lq
        #nh = len([i for i in diffs if i >= hq]) # high
        #nl = len([i for i in diffs if i <= lq]) # low
        nh = len([i for i in diff_array if i >= hq]) # high
        nl = len([i for i in diff_array if i <= lq]) # low
        # store values
        stats['hq'] = hq
        stats['lq'] = lq
        stats['nh'] = nh
        stats['nl'] = nl
        #stats['diffs'] = diffs
        stats['diffs'] = diff_array
        # store stats by sample
        data[sample] = stats
        genediffs[sample] = diff_dict

    return data, genediffs


def build_plots(data, sigval, target_sample, outdir):
    
    # define colors
    bcolor = '#969696'
    lcolor = '#000000'
    mcolor = '#525252'
    gridM = '#bdbdbd'
    alpha = 0.3

    for sample, stats in data.items():

        print(f'\t\tPlotting Sample: {sample} ...')

        # define outfile
        outfile = f'{outdir}/{target_sample}-{sample}_diffcov_plot.png'

        # setup to plot
        fig, ax = plt.subplots(figsize=(8.5, 4), dpi=300)
        # plot histogram
        ax.hist(stats['diffs'], bins=30, color=bcolor, alpha=alpha)

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

        # plot sigval line
        axes = plt.gca()
        y_min, y_max = axes.get_ylim()
        ytextpos = y_max - (y_max / 20)

        ax.axvline(
            x=stats['lq'], ymin=0, ymax=1, color=mcolor, linewidth=2,
            linestyle='--', label='Lower Quantile'
            )
        ax.axvline(
            x=stats['hq'], ymin=0, ymax=1, color=mcolor, linewidth=2,
            linestyle='--', label='Upper Quantile'
            )
        qtext = (
            f" Sig. alpha:  {sigval}\n"
            f" Sig. High DiffCov: {stats['hq']}\n"
            f" Sig. Low DiffCov: {stats['lq']}\n"
            f" Sig. High Genes: {stats['nh']}\n"
            f" Sig. Low Genes: {stats['nl']}"
            )
        ax.text(
            stats['lq'], ytextpos, qtext, color=mcolor, fontsize=11,
            horizontalalignment='left', verticalalignment='top'
            )

        # set plot title and axis labels
        ax.set_title(
            f'{target_sample}-{sample}',
            fontsize=18, fontweight='bold', y=1.02
            )
        ax.set_xlabel(
            f'log(TAD/AVG) Difference',
            fontsize=14, fontweight='bold', y=-0.02
            )
        ax.set_ylabel(
            'Gene Count',
            fontsize=14, fontweight='bold', x=-0.02
            )

        # save figure and close
        plt.tight_layout()
        plt.savefig(outfile)
        plt.close()

    return True


def write_output(anno, fasta, diffs, data, outdir, target_sample, nmtd, td):

    # iterage through samples
    for sample, genediffs in diffs.items():

        # setup outfiles
        base = f'{outdir}/{target_sample}-{sample}'
        anno_high = open(f'{base}_diffcov_high.tsv', 'w')
        anno_low = open(f'{base}_diffcov_low.tsv', 'w')

        if anno:
            outheader = (
                f"Gene\tlog(TAD/AVG) Diff\tTarget_Norm_TAD\tTest_Norm_TAD\t"
                f"Target_TAD\tTest_TAD\t{anno['Header']}\n"
                )
        else:
            outheader = (
                "Gene\tlog(TAD/AVG) Diff\tTarget_Norm_TAD\t"
                "Test_Norm_TAD\tTarget_TAD\tTest_TAD\n"
                )

        anno_high.write(outheader)
        anno_low.write(outheader)

        if fasta:
            fasta_high = open(f'{base}_diffcov_high.fa', 'w')
            fasta_low = open(f'{base}_diffcov_low.fa', 'w')

        # set high/low values
        hq = data[sample]['hq']
        lq = data[sample]['lq']

        # iterate genes >= hq and write files.
        high_genes = {k: v for k, v in genediffs.items() if v >= hq}
        for gene, tad in high_genes.items():
            if fasta: fasta_high.write(f'>{gene}\n{fasta[gene]}\n')
            nmtdA = nmtd[target_sample][gene] # target sample value
            tdA = td[target_sample][gene] # target sample value
            nmtdB = nmtd[sample][gene] #test sample value
            tdB = td[sample][gene] #test sample value
            if anno:
                hout = f'{gene}\t{tad}\t{nmtdA}\t{nmtdB}\t{tdA}\t{tdB}\t{anno[gene]}\n'
                anno_high.write(hout)
            else:
                hout = f'{gene}\t{tad}\t{nmtdA}\t{nmtdB}\t{tdA}\t{tdB}\n'
                outline = anno_high.write(hout)

        # iterate genes <= lq and write files.
        low_genes = {k: v for k, v in genediffs.items() if v <= lq}
        for gene, tad in low_genes.items():
            if fasta: fasta_low.write(f'>{gene}\n{fasta[gene]}\n')
            nmtdA = nmtd[target_sample][gene] # target sample value
            tdA = td[target_sample][gene] # target sample value
            nmtdB = nmtd[sample][gene] #test sample value
            tdB = td[sample][gene] #test sample value
            if anno:
                lout = f'{gene}\t{tad}\t{nmtdA}\t{nmtdB}\t{tdA}\t{tdB}\t{anno[gene]}\n'
                anno_low.write(lout)
            else:
                lout = f'{gene}\t{tad}\t{nmtdA}\t{nmtdB}\t{tdA}\t{tdB}\n'
                outline = anno_low.write(lout)

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_directory',
        help='Please specify Coverage Magic output directory!',
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
        '-o', '--output_directory',
        help='Please specify the output directory!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-s', '--significance_alpha',
        help='OPTIONAL: Change significance alhpa (Default=0.025)!',
        metavar='',
        type=float,
        default=0.025,
        required=False
        )
    parser.add_argument(
        '-a', '--annotation_directory',
        help='OPTIONAL: Specify the annotations file from MicrobeAnnotator!',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-f', '--fasta_directory',
        help='OPTIONAL: Specify the genes fasta file!',
        metavar='',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script ...\n')

    # define input parameters
    indir = args['input_directory']
    target_sample = args['target_sample']
    outpre = args['output_directory']
    sigval = args['significance_alpha']
    annotations = args['annotation_directory']
    genefastas = args['fasta_directory']

    # If annotations parse annotations files
    if annotations:
        annotations = parse_annotations(annotations)
        print('\n\nFinished parsing annotations. Continuing ...\n')

    # If fasta parse fasta files
    if genefastas:
        genefastas = parse_fasta(genefastas)
        print('\n\nFinished parsing genes fasta. Continuing ...\n')

    # Parse the input file
    lgnmtd, nmtd, td = parse_input(indir)
    print('\n\nFinished parsing input files. Analyzing data ...\n')

    # Analyze data
    data, genediffs = analyze_data(lgnmtd, sigval, target_sample)
    print('\n\nFinished analyzing data. Building plots ...\n')

    # Build and write out the plots.
    _ = build_plots(data, sigval, target_sample, outpre)
    print('\n\nFinished building plots. Writing output files ...')

    # Write output files
    _ = write_output(
            annotations, genefastas, genediffs,
            data, outpre, target_sample,
            nmtd, td
            )

    print('\n\nComplete success space cadet!! We finished without errors.\n\n')


if __name__ == "__main__":
    main()
