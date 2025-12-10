#!/usr/bin/env python

''' Calculate Percent reads mapped to MAGs

Takes two files:
    1) Coupled Metagenome
    2) Concatenated Tabular Blast

Returns the number of reads mapped in blast file above 90% pident
(default) divided by the total number of reads in the metagenome file.

i.e. The percent of the community represented by the MAG collection.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: July 13 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
from pathlib import Path

def meta_len(file):
    with open(file) as f:
        for i, l in enumerate(f, 1):
            pass
    return i

def blast_len(file, cutoff):
    i = 0
    with open(file) as f:
        for l in f:
            X = l.split('\t')
            pid = float(X[2])
            if pid >= cutoff: i += 1
    return i

def percent_community_represented(blastfile, metafile, cutoff):

    # number of reads in metagenome fasta
    # count number of lines and divide by 2
    metalen = meta_len(metafile) / 2

    # number of reads passing cutoff in blast file
    blastlen = blast_len(blastfile, cutoff)

    community_represented = blastlen / metalen * 100

    metaname = metafile.split('/')[-1].split('.')[0]
    print(f'{metaname}\t{community_represented:.2f}')

def main():

        # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-b', '--filtered_blast_file',
        help='Please specify the tabular blast file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--corresponding_metagenome_file',
        help='Please specify the corresponding metagenome file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--pident_cutoff',
        help='pident value for population cutoff! (default = 90.0) ',
        metavar=':',
        type=float,
        required=False,
        default=90.0
        )
    args=vars(parser.parse_args())

    percent_community_represented(
                    args['filtered_blast_file'],
                    args['corresponding_metagenome_file'],
                    args['pident_cutoff']
                    )

if __name__ == "__main__":
    main()