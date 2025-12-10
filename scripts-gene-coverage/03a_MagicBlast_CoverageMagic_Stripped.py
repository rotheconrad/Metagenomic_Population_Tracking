#!/usr/bin/env python

'''Calculates ANIr and Sequence Coverage from Magic Blast Output.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Magic Blast output should be filtered prior to using this script   !!
!! Use 01c_ShortRead_Filter.py or other method.                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This script calculates ANIr and sequence coverage (as depth and    !!
!! breadth) from tabular Magic Blast output for the whole genome or   !!
!! or MAG, each contig in the genomic sequence fasta file, each       !!
!! predicted coding sequence (protein coding gene), and each          !!
!! intergenic region. It requires the metagenomic fasta file used as  !!
!! the blast queries, the genonimic reference sequence used as the    !!
!! database or subject, and a fasta file of predicted protein coding  !!
!! sequences predicted by Prodigal from the genomic reference.        !!
!! The *_CDS_from_genomic.fna files from NCBI are usually in order.   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Coverage calculated as Truncated Average Depth (TAD):
    * Set TAD to 100 for no truncatation.
    * TAD 80 removes the top 10% and bottom 10% of base pair depths and
      caluclates coverage from the middle 80% of values. Intended to 
      reduce effects of conserved motif peaks and contig edge valleys.
    * Coverage = base pairs recruited / length of genome, contig, or gene

Coverage calculated as Breadth:
    * number of positions in reference sequence covered by at least
      one read alignment divided the length of the reference sequence.

Relative Abundance is calculated as:
    * base pairs recruited / base pairs in metagenome * 100
    * It is the percent of base pairs recruited out of the total
      base pairs sequenced in the metagenome.

ANIr is calculated as:
    * average percent identity of sequence alignments for all reads 
      (should be 1 blast match per read)

This tool takes the following input parameters:

    * Tabular Blast file containing results for 1 genome and 1 metagenome
    * *_genomic.fna file from NCBI used as reference for blast search.
    * Metagenome fastq (or fasta) file used as queries for blast search.
    * *_CDS_from_genomic.fna file of predicted genes for *_genomic.fna.

This script returns the following files:

    * 3 column tsv output of Contig(or gene_name), coverage(or ANIr), 
      sequence length.
    * Writes 8 files total:
        - {out_file_prefix}_genome.tsv
        - {out_file_prefix}_contig_tad.tsv
        - {out_file_prefix}_contig_breadth.tsv
        - {out_file_prefix}_contig_anir.tsv
        - {out_file_prefix}_gene_tad.tsv
        - {out_file_prefix}_gene_breadth.tsv
        - {out_file_prefix}_intergene_tad.tsv
        - {out_file_prefix}_intergene_breadth.tsv

*_gene_* files contain values for the CDS regions.
*_intergene_* files contain values for the inter-CDS regions.

This script requires the following packages:

    * argparse, sys
    * collection.defaultdict
    * itertools

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February 12th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, sys
import itertools
from collections import defaultdict


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


def read_genome_lengths(rgf):
    """ Reads genome lengths file returns dict genome_name: length """

    rgf_tad = defaultdict(dict) # initialize dicts
    rgf_len = {}
    wg_len = 0

    # read through genome fasta and build dictionary of dictionary
    # Containing base pair position for length of each contig.
    with open(rgf, 'r') as f:
        for name, seq in read_fasta(f):

            contig_name = name.split(' ')[0][1:]
            length = len(seq) # calculate length of contig

            rgf_len[contig_name] = length
            wg_len += length

            # This populates the dictionary with value of zero for each
            # base pair position for each contig in the genome fasta
            for i in range(1, length+1, 1):
                rgf_tad[contig_name][i] = 0

    return rgf_tad, rgf_len, wg_len


def calc_genome_coverage(tbf, rgf_tad, lthd, uthd):
    """ Reads tabblast file and adds coverage by genome position """

    read_count = 0
    bpmap = 0

    with open(tbf, 'r') as f:
        for l in f:
            # Progress Tracker
            if read_count % 50000 == 0:
                print(f'... Blast matches processed so far ... {read_count:013}')
            read_count += 1

            # Skip magic blast header
            if l.startswith('#'): continue

            # split each line and define variables of interest
            X = l.rstrip().split('\t')
            pident = float(X[2])
            contig_name = X[1]
            strt = min(int(X[8]), int(X[9]))
            stp = max(int(X[8]), int(X[9]))
            qlen = int(X[15])
            bpmap += qlen
            
            # for each read above the user specified threshold, add
            # coverage of +1 to each basepair position for length of
            # read alignment along the subject sequence.

            if pident > lthd and pident < uthd:

                # add +1 coverage to each position the read covers
                for i in range(strt, stp+1, 1):
                    rgf_tad[contig_name][i] += 1

    print(f'... Total Blast matches process: {read_count:013}')

    return rgf_tad, bpmap


def truncate(x, tad):
    """ returns tad range of a list/array """

    xsorted = sorted(x)
    xlen = len(x) # get length of the list
    inverse_tad = round((1.0 - tad) / 2.0, 2) # to get top and bottom 
    q = int(xlen * inverse_tad) # to get top and bottom
    bottom = q
    top = xlen - q
    tad_range = xsorted[bottom:top] # slice list

    #print(xlen, bottom, top, len(tad_range), sum(xsorted[:bottom+1]))

    return tad_range


def get_contig_tad(rgf_tad, tad):
    """ reads through rgf_tad and returns dict of tads by contig """

    contig_tad = {}
    contig_breadth = {}
    wg_tad = []

    for k, v in rgf_tad.items():
        values = list(v.values())
        breadth = sum(i > 0 for i in values) / len(values)
        contig_breadth[k] = breadth
        coverage = get_average(values, tad)
        contig_tad[k] = coverage
        wg_tad.extend(values)

    return contig_tad, contig_breadth, wg_tad


def get_gene_tad(gn_tad, tad):
    """ reads through gn_tad and returns dict of tads by gene """

    gene_tad = {}
    gene_breadth = {}

    for k, v in gn_tad.items():
        breadth = sum(i > 0 for i in v) / len(v)
        gene_breadth[k] = breadth
        gene_tad[k] = get_average(v, tad)

    return gene_tad, gene_breadth
    

def get_average(l, tad):
    """ returns anir from truncated list """

    trunc_val = truncate(l, tad)

    if sum(trunc_val) == 0:
        average = 0
    else:
        average = sum(trunc_val) / len(trunc_val)

    return average


def write_file(tad, breadth, length, outpre, outpost, precision):
    """ writes dictionary to file """
    
    outfile = outpre + outpost

    with open(outfile, 'w') as o:
        o.write('Contig_Name\tTAD80\tLength\tBreadth\n')
        for k, v in tad.items():
            o.write(
                f'{k}\t{v:.{precision}f}\t{length[k]}\t'
                f'{breadth[k]:.{precision}f}\n'
                )

    return True


def calc_contig_stats(
    rgf_tad, rgf_len, tad, outpre, precision
    ):
    """Calculate ANIr, TAD and breadth and write to files for Contigs"""

    print('... Calculating TADs for Contigs.')
    contig_tad, contig_breadth, wg_tad = get_contig_tad(rgf_tad, tad)

    print('... Writing Contig file.')
    _ = write_file(
                contig_tad,
                contig_breadth,
                rgf_len,
                outpre,
                '_contig.tsv',
                precision
                )

    contig_tad = None
    contig_breadth = None

    return wg_tad


def operator(rgf, tbf, lthd, uthd, tad, outpre, precision):
    """ Runs the different functions and writes out results """

    tadp = tad / 100

    print(
        f'Using values {tad}% for TAD & {lthd}%, {uthd}% for pIdent Threshold'
        )

    print('\nPreparing base pair array for each contig in genome.')
    rgf_tad, rgf_len, wglen = read_genome_lengths(rgf)

    print(
        '\nCalculating coverage for each base pair position in genome. \n'
        'This can take a while depending on the number of blast results.'
        )
    rgf_tad, bpmap = calc_genome_coverage(tbf, rgf_tad, lthd, uthd)

    print(
        f'\nCalculating {tad}% truncated average depth '
        f'with {lthd}% lower and {uthd}% upper percent sequence identity.'
        )

    wg_tad = calc_contig_stats(
                    rgf_tad, rgf_len, tadp, outpre, precision
                    )

    # Clear from memory
    rgf_tad = None
    rgf_len = None

    print('\nScript seems to have finished successfully.\n')


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-g', '--ref_genome_file',
        help='Please specify the genome fasta file!',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-b', '--tabular_blast_file',
        help='Please specify the tabular blast file!',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--precision',
        help=
            '(Optional) Change decimal places in output (default = 2).',
        #metavar='',
        type=int,
        required=False,
        default=2
        )
    parser.add_argument(
        '-c', '--pIdent_threshold_lower',
        help=
            '(Optional) Lower percent sequence identity of reads to include '
            'coverage calculations (Default = 94.99). [pIdent > value].',
        #metavar='',
        type=float,
        required=False,
        default=94.99
        )
    parser.add_argument(
        '-u', '--pIdent_threshold_upper',
        help=
            '(Optional) Upper percent sequence identity of reads to include '
            'coverage calculations (Default = 100.01). [pIdent < value].',
        #metavar='',
        type=float,
        required=False,
        default=100.01
        )
    parser.add_argument(
        '-d', '--truncated_avg_depth_value',
        help='(Optional) Specify a different TAD value! (Default = 80)',
        #metavar='',
        type=float,
        required=False,
        default=80
        )
    parser.add_argument(
        '-o', '--out_file_prefix',
        help='What do you like the output file prefix to be?',
        #metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')
    operator(
            args['ref_genome_file'],
            args['tabular_blast_file'],
            args['pIdent_threshold_lower'],
            args['pIdent_threshold_upper'],
            args['truncated_avg_depth_value'],
            args['out_file_prefix'],
            args['precision'],
            )


if __name__ == "__main__":
    main()
