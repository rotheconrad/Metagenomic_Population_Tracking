#!/usr/bin/env python

'''Finds genes in the lowest 2.5% range of gene coverage.

This scripts takes the *_gene_tad.tsv output from the coverage magic
script and finds genes with coverage below the 0.025 quantile.

This equates to statistically significant genes using an alpha of 0.025

It plots a histogram with the 2.5 percentile line and returns a list
of genes below this value.

Quantile and percentile represent the same thing
e.g. 0.025 quantile = 2.5% percentile

Gene coverage values are not necessarily normally distributed and so
the non-parametric quantile/percentile metric is used to determine which
genes have significantly lower coverage than them majority of genes.

Input:
    - *_gene_tad.tsv from coverage magic script.
    - genes.faa fasta file of the genes (optional).
    - *.annotations file from MicrobeAnnotator (optional).

Output:
    - Histogram of gene coverage values as .png
    - tsv file of significantly low coverage genes.
    - fasta file with significantly low coverage genes (optional).

Optional: If an accompanying microbe annotator *.annotations file (found
in the annotation_results output folder) is provided this scripts
significantly low coverage gene output will include the matching
annotation.

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


def parse_gene_tad(infile):
    
    genetads = {}

    with open(infile, 'r') as file:
        header = file.readline()

        for line in file:
            X = line.rstrip().split('\t')
            gene = X[0]
            tad = float(X[1])
            genetads[gene] = tad

    return genetads


def find_significant_low(genetads, sigval):
     
     data = {}

     # create array of measured tad values
     tads = [v for k, v in genetads.items()]

     # find the value of the sigval percentile
     q = np.quantile(tads, sigval)

     # how many genes below this value?
     n = len([i for i in tads if i <= q])

     # how many genes == 0?
     z = len([i for i in tads if i == 0])

     data['tads'] = tads
     data['q'] = q
     data['n'] = n
     data['z'] = z

     return data


def plot_histogram(data, sigval, outpre):
    
    # define colors
    #acolor = '#252525'
    bcolor = '#969696'
    lcolor = '#000000'
    mcolor = '#525252'
    gridM = '#bdbdbd'
    alpha = 0.3

    # setup to plot
    fig, ax = plt.subplots(figsize=(8.5, 4), dpi=300)
    # plot histogram
    ax.hist(data['tads'], bins=30, color=bcolor, alpha=alpha)

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
        x=data['q'], ymin=0, ymax=1, color=mcolor, linewidth=2,
        linestyle='--', label='Lower Quantile'
        )
    qtext = (
        f" Sig. alpha:  {sigval}\n"
        f" Sig. Coverage: {data['q']}\n"
        f" Genes Below: {data['n']}\n"
        f" Genes Zero: {data['z']}"
        )
    ax.text(
        data['q'], ytextpos, qtext, color=mcolor, fontsize=11,
        horizontalalignment='left', verticalalignment='top'
        )

    # set plot title and axis labels
    ax.set_title(
        f'Distribution of Gene Coverage Depth',
        fontsize=18, fontweight='bold', y=1.02
        )
    ax.set_xlabel(
        f'Coverage Depth',
        fontsize=14, fontweight='bold', y=-0.02
        )
    ax.set_ylabel(
        'Gene Count',
        fontsize=14, fontweight='bold', x=-0.02
        )

    # save figure and close
    plt.tight_layout()
    plt.savefig(f'{outpre}_gene_coverage.png')
    plt.close()

    return True


def write_output(annotations, genefasta, genetads, data, outpre):
    
    q = data['q']

    sig_genes = {k: v for k, v in genetads.items() if v <= q}

    if genefasta:
        gfasta = open(f'{outpre}.fasta', 'w')

    outfile = open(f'{outpre}_sig_genes.tsv', 'w')

    if annotations:
        outheader = f"Gene\tTAD\t{annotations['Header']}\n"
    else:
        outheader = "Gene\tTAD\n"

    outfile.write(outheader)

    for gene, tad in sig_genes.items():

        if genefasta:
            lineout = f'>{gene}\n{genefasta[gene]}\n'
            gfasta.write(lineout)

        if annotations:
            outline = f'{gene}\t{tad}\t{annotations[gene]}\n'
        else:
            outline = f'{gene}\t{tad}\t'

        outfile.write(outline)

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_gene_tad_file',
        help='Please specify the *_gene_tad.tsv input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the filename prefix for the output files!',
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
        '-a', '--annotation_file',
        help='OPTIONAL: Specify the annotations file from MicrobeAnnotator!',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-f', '--fasta_file',
        help='OPTIONAL: Specify the genes fasta file!',
        metavar='',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script ...\n')

    # define input parameters
    infile = args['input_gene_tad_file']
    outpre = args['output_prefix']
    sigval = args['significance_alpha']
    annotations = args['annotation_file']
    genefasta = args['fasta_file']

    # If annotations parse annotations file
    if annotations:
        annotations = parse_annotations(annotations)
        print('\n\nFinished parsing annotations. Continuing ...\n')

    # If fasta parse fasta file
    if genefasta:
        genefasta = parse_fasta(genefasta)
        print('\n\nFinished parsing genes fasta. Continuing ...\n')

    genetads = parse_gene_tad(infile)
    print('\n\nFinished parsing gene tad file. Continuing ...\n')

    data = find_significant_low(genetads, sigval)
    print('\n\nFinished low gene significance test. Plotting ...\n')

    _ = plot_histogram(data, sigval, outpre)
    print('\n\nFinished plotting. Writing output files ...\n')

    _ = write_output(annotations, genefasta, genetads, data, outpre)
    print('\n\nComplete success space cadet!! We finished without errors.\n\n')


if __name__ == "__main__":
    main()
