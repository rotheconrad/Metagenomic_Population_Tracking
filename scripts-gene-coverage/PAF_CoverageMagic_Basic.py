#!/usr/bin/env python

'''Calculates ANIr and Sequence Coverage from PAF format file.

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

    Required:
        * PAF file containing results for 1 genome and 1 metagenome.
        * MAG or genome fasta file the reads were mapped to.
    Optional:
        * Metagenome fasta file to compute relative abundance. 
        * Prodigal CDS fasta file to compute gene and intergene TAD.

This script returns the following files:

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


def calc_genome_coverage(paf_file, rgf_tad, min_pid, min_frac, min_qual):
    """ Reads and filters PAF file and adds coverage by genome position """

    rgf_anir = defaultdict(list)
    read_count = 0

    data = {}
    total_match = 0
    dup_match = 0
    passing_match = 0

    with open(paf_file, 'r') as file:
        for line in file:
            if line[0] == '@': continue # skip the header
            total_match += 1
            X = line.split(sep='\t')
            qname = X[0].split(' ')[0] # query sequence name
            contig_name = X[5].split(' ')[0] # target sequence name
            if qname == contig_name: continue # remove any self matches
            strt = min(int(X[7]), int(X[8])) # target start
            stp = max(int(X[7]), int(X[8])) # target end
            qlen = int(X[1]) # length of the query sequence (read length)
            matches = int(X[9]) # number of matches across alignment length
            alen = int(X[10]) # length of the alignment
            pid = 100*(matches / alen)
            qual = int(X[11]) # mapping quality
            afrac = alen / qlen # aligned fractions
            # some quality and best match filtering
            # we aren't interesting in anything below 70% pid generally
            if afrac >= min_frac and pid >= min_pid and qual >= min_qual:
                # count each query sequence only once.
                # For duplicate mathces keep the one with the highest qual
                if qname in data:
                    dup_match += 1
                    old_qual = data[qname][0]
                    if qual > old_qual:
                        data[qname] = [qual, pid, strt, stp, contig_name]
                else:
                    passing_match += 1
                    data[qname] = [qual, pid, strt, stp, contig_name]

    print(
        f'\n\t\tFile Name: {paf_file}\n'
        f'\t\tTotal matches in file: {total_match}\n'
        f'\t\tDuplicate matches: {dup_match}\n'
        f'\t\tTotal passing matches: {passing_match}\n'
        )
            
    # for each read above the user specified threshold, add
    # coverage of +1 to each basepair position for length of
    # read alignment along the subject sequence.

    for read, metrics in data.items():
        pid = metrics[1]
        strt = metrics[2]
        stp = metrics[3]
        contig_name = metrics[4]

        rgf_anir[contig_name].append(pid) # add pident to anir list
        # add +1 coverage to each position the read covers
        for i in range(strt+1, stp+2, 1):
            rgf_tad[contig_name][i] += 1

    return rgf_tad, rgf_anir


def retrieve_prodigal_gene_coverage(pgf, rgf_tad):
    """ Retrieves list of depths for each bp position of gene length """

    gn_tad = defaultdict(list) # initialize dictionary
    gn_len = {}

    intergn_tad = defaultdict(list) # initialize dictionary
    intergn_len = {}

    with open(pgf, 'r') as f:
        stp = 1 # initial stp value for intergene calculation
        for name, seq in read_fasta(f):
            X = name.split(' # ')
            gene_name = X[0][1:]
            contig_name = '_'.join(gene_name.split('_')[:-1])

            strt = min(int(X[1]), int(X[2]))

            # Define intergenic or between CDS regions
            intergene_strt = stp + 1 # start of inter-CDS region
            intergene_stp = strt # stop of inter-CDS region
            intergene_len = intergene_stp - intergene_strt

            stp = max(int(X[1]), int(X[2]))

            intergene_name = (
                f'{contig_name}_intergene_{intergene_strt}-{intergene_stp}'
                )

            gn_len[gene_name] = len(seq)
            intergn_len[intergene_name] = intergene_len

            # Get depth values for gene (CDS) regions
            for i in range(strt, stp+1, 1):
                gn_tad[gene_name].append(rgf_tad[contig_name][i])

            # Get depth values for intergene (inter-CDS) regions
            for i in range(intergene_strt, intergene_stp+1, 1):
                intergn_tad[intergene_name].append(rgf_tad[contig_name][i])

        # Get intergene region after last predicted coding region.
        intergene_strt = stp + 1
        intergene_stp = len(rgf_tad[contig_name])
        intergene_len = intergene_stp - intergene_strt

        intergene_name = (
            f'{contig_name}_intergene_{intergene_strt}-{intergene_stp}'
            )
        intergn_len[intergene_name] = intergene_len
        # Get depth values for intergene (inter-CDS) regions
        for i in range(intergene_strt, intergene_stp+1, 1):
            intergn_tad[intergene_name].append(rgf_tad[contig_name][i])

    return gn_tad, gn_len, intergn_tad, intergn_len


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
    

def get_contig_anir(rgf_anir, tad):
    """ loops through in_d and calculates tad/ani for each key """

    contig_anir = {}
    wg_anir = []

    for k, v in rgf_anir.items():
        average = get_average(v, tad)
        #if average > 0: contig_anir[k] = average
        contig_anir[k] = average
        wg_anir.extend(v)

    return contig_anir, wg_anir


def get_average(l, tad):
    """ returns anir from truncated list """

    trunc_val = truncate(l, tad)

    if sum(trunc_val) == 0:
        average = 0
    else:
        average = sum(trunc_val) / len(trunc_val)

    return average


def get_relative_abundance(wg_tad, mtg):
    """ calculates and returns relative abundance from wg_TAD """

    total_metagenome_bp = 0

    # check if fasta or fastq
    file_type = mtg.split('.')[-1]
    fqtype = ['fastq', 'fq']
    fatype = ['fasta', 'fna', 'fst', 'fa']

    if file_type in fqtype:
        line_count = 0

        with open(mtg, 'r') as f:
            for l in f:
                line_count += 1
                if line_count%4 == 0:
                    total_metagenome_bp += len(l.rstrip())
    elif file_type in fatype:
        with open(mtg, 'r') as f:
            for name, seq in read_fasta(f):
                total_metagenome_bp += len(seq)
    else:
        print(
            'Error in determining metagenome format of fasta or fastq. '
            'Please double check metagenome file type and try again. '
            'Metagenome file should be either fasta or fastq format with '
            f'file extension of one of {fqtype} or {fatype}.\n\n'
            'If there is a file extension you would like added, please '
            'submit a feature request through the issues tab of the '
            'GitHub repo at: '
            'https://github.com/rotheconrad/00_in-situ_GeneCoverage/issues'
            '.\n I will be happy to add aditional file extensions.\n\n'
            )
        sys.exit()

    relabndc = (sum(wg_tad) / total_metagenome_bp) * 100

    return relabndc, total_metagenome_bp


def write_file(in_d, len_d, outpre, outpost, precision):
    """ writes dictionary to file """
    
    outfile = outpre + outpost

    with open(outfile, 'w') as o:
        o.write('Name\tValue\tLength\n')
        for k, v in in_d.items():
            o.write(f'{k}\t{v:.{precision}f}\t{len_d[k]}\n')


def calc_contig_stats(
    rgf_tad, rgf_anir, rgf_len, tad, outpre, precision
    ):
    """Calculate ANIr, TAD and breadth and write to files for Contigs"""

    print('... Calculating TADs for Contigs.')
    contig_tad, contig_breadth, wg_tad = get_contig_tad(rgf_tad, tad)

    print('... Writing Contig TAD file.')
    _ = write_file(contig_tad, rgf_len, outpre, '_contig_tad.tsv', precision)

    print('... Writing Contig Breadth file.')
    _ = write_file(
        contig_breadth, rgf_len, outpre, '_contig_breadth.tsv', precision
        )

    contig_tad = None
    contig_breadth = None

    print('... Calculating ANIr for Contigs')
    contig_anir, wg_anir = get_contig_anir(rgf_anir, tad)

    print('... Writing Contig ANIr file.')
    _ = write_file(contig_anir, rgf_len, outpre, '_contig_anir.tsv', precision)

    contig_anir = None

    return wg_tad, wg_anir


def calc_gene_stats(
    gn_tad, gn_len, tad, outpre, precision
    ):
    """Calculate ANIr, TAD and breadth and write to files for Genes"""

    print('... Calculating TADs for Genes')
    gene_tad, gene_breadth = get_gene_tad(gn_tad, tad)

    print('... Writing Gene TAD file.')
    _ = write_file(gene_tad, gn_len, outpre, '_gene_tad.tsv', precision)

    print('... Writing Gene Breadth file.')
    _ = write_file(gene_breadth, gn_len, outpre, '_gene_breadth.tsv', precision)

    gene_tad = None
    gene_breadth = None


def calc_intergene_stats(
    intergn_tad, intergn_len, tad, outpre, precision
    ):

    """Calculate ANIr, TAD and breadth and write to Intergene files"""

    print('... Calculating TADs for Inter-Gene Regions')
    intergene_tad, intergene_breadth = get_gene_tad(intergn_tad, tad)

    print('... Writing Intergene TAD file.')
    _ = write_file(
        intergene_tad, intergn_len, outpre, '_inter-gene_tad.tsv', precision
        )

    print('... Writing Intergene Breadth file.')
    _ = write_file(
        intergene_breadth, intergn_len, outpre,
        '_inter-gene_breadth.tsv', precision
        )

    intergene_tad = None
    intergene_breadth = None


def calc_genome_stats(
    mtg, wg_tad, wg_anir, wglen, tad, min_pid, outpre, precision
    ):
    """Calculates ANIr, TAD and breadth and writes to Genome files"""

    print('... Calculating TAD for Genome')
    wgbreadth = sum(i > 0 for i in wg_tad) / len(wg_tad)
    wgtad = get_average(wg_tad, tad)

    if mtg is None:
        relabndc = 'n/a'
        total_metagenome_bp = 'n/a'
    elif mtg.isdigit():
        total_metagenome_bp = int(mtg)
        relabndc = (sum(wg_tad) / total_metagenome_bp) * 100
    else:
        print('... Calculating Total Metagenome Size & Relative Abundance')
        relabndc, total_metagenome_bp = get_relative_abundance(wg_tad, mtg)
        relabndc = f'{relabndc:.{precision}f}%'

    print('... Calculating ANI for Genome')
    wganir = get_average(wg_anir, tad)

    wg_header = (
            f'Genome_Name\tTAD_{tad*100}\tBreadth\tANIr_{min_pid}\t'
            f'Relative_Abundance(%)\tGenome_Length(bp)\tMetagenome_Length(bp)\n'
            )
    wg_lineout = (
            f'{outpre}\t{wgtad:.{precision}f}\t{wgbreadth:.{precision}f}\t'
            f'{wganir:.{precision}f}%\t'
            f'{relabndc}\t{wglen}\t{total_metagenome_bp}\n'
            )

    with open(f'{outpre}_genome.tsv', 'w') as o:
        o.write(wg_header)
        o.write(wg_lineout)

    print('\nWhole Genome Values:\n')
    print(wg_header[:-1])
    print(wg_lineout[:-1], '\n')


def operator(
    ref_genome, paf_file, gene_CDS, metaG, min_pid, min_frac, min_qual, tad, outpre
    ):
    """ Runs the different functions and writes out results """

    tadp = tad / 100
    precision = 2 # number of decimals places to keep.

    print(
        f'Using {min_pid}% sequence identity as the read cutoff\n'
        f'Computing {tad}% truncated average sequencing depth (TAD)'
        )

    print('\nPreparing base pair array for each contig in genome.')
    rgf_tad, rgf_len, wglen = read_genome_lengths(ref_genome)

    print(
        '\nCalculating coverage for each base pair position in genome. \n'
        'This can take a while depending on the number of blast results.'
        )
    rgf_tad, rgf_anir = calc_genome_coverage(
                                paf_file, rgf_tad, min_pid, min_frac, min_qual
                                )

    ### Check for Prodigal or NCBI. #####################################
    if gene_CDS:
        print(
            '\nRetrieving coverage for each contig & gene.'
            )

        (
        gn_tad, gn_len, intergn_tad, intergn_len
            ) = retrieve_prodigal_gene_coverage(gene_CDS, rgf_tad)

        do_genes = True

    else:
        print(
            '\n\n!! No gene prediction file entered.\n'
            '!! NOT Calculating values by genes.\n'
            '!! Calculating values for by contig and whole genome only.\n'
            )

        do_genes = False

    ### End Predicted Gene File Check ###################################

    print(
        f'\nCalculating {tad}% truncated average depth '
        f'and {min_pid}% threshold ANIr'
        )

    wg_tad, wg_anir = calc_contig_stats(
                    rgf_tad, rgf_anir, rgf_len, tadp, outpre, precision
                    )

    # Clear from memory
    rgf_tad = None
    rgf_anir = None
    rgf_len = None

    if do_genes == True:

        _ = calc_gene_stats(
            gn_tad, gn_len, tadp, outpre, precision
            )

        # Clear from memory
        gn_tad = None
        gn_len = None

        _ = calc_intergene_stats(
            intergn_tad, intergn_len, tadp, outpre, precision
            )

        intergn_tad = None
        intergn_len = None

    _ = calc_genome_stats(
        metaG, wg_tad, wg_anir, wglen, tadp, min_pid, outpre, precision
        )

    print('\nScript seems to have finished successfully.\n')


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-ref', '--ref_genome_file',
        help='Please specify the reference genome fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-paf', '--paf_format_file',
        help='Please specify the PAF format file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-cds', '--prodigal_gene_fasta',
        help=
            '(Optional) Use this option to report values for genes and '
            'intergenic regions predicted with Prodigal.',
        metavar='',
        type=str,
        required=False,
        default=None
        )
    parser.add_argument(
        '-mtg', '--metagenome_fasta',
        help=
            '(Optional) To calculate relative abundance specify either the path'
            ' to the metagenome file, or the size of the metagenome in base '
            'pairs.',
        metavar='',
        type=str,
        required=False,
        default=None
        )
    parser.add_argument(
        '-min_pid', '--minimum_sequence_identity',
        help='(OPTIONAL) Minimum percent sequence identity (Default = 95).',
        metavar='',
        type=int,
        required=False,
        default=95
        )
    parser.add_argument(
        '-min_frac', '--minimum_aligned_frac',
        help='(OPTIONAL) Specify the minimum aligned_frac (Default = 0.5).',
        metavar='',
        type=float,
        required=False,
        default=0.5
        )
    parser.add_argument(
        '-min_qual', '--minimum_quality',
        help='(OPTIONAL) Specify the minimum PAF quality score (Default = 5).',
        metavar='',
        type=int,
        required=False,
        default=5
        )
    parser.add_argument(
        '-tad', '--truncated_avg_depth',
        help='(Optional) Specify a different TAD value! (Default = 80)',
        metavar='',
        type=float,
        required=False,
        default=80
        )
    parser.add_argument(
        '-out', '--out_file_prefix',
        help='What do you like the output file prefix to be?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # define input params
    ref_genome = args['ref_genome_file']
    paf_file = args['paf_format_file']
    gene_CDS = args['prodigal_gene_fasta']
    metaG = args['metagenome_fasta']
    min_pid = args['minimum_sequence_identity']
    min_frac = args['minimum_aligned_frac']
    min_qual = args['minimum_quality']
    tad = args['truncated_avg_depth']
    outpre = args['out_file_prefix']

    # Do what you came here to do:
    print('\n\nRunning Script...\n')
    operator(
            ref_genome,
            paf_file,
            gene_CDS,
            metaG,
            min_pid,
            min_frac,
            min_qual,
            tad,
            outpre
                )


if __name__ == "__main__":
    main()
