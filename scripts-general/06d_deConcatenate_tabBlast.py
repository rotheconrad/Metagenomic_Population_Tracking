#!/usr/bin/env python

''' De-concatenate Tabular Blast File

Parses unique sample ID in second column and writes new output files
for each unique ID found. Second column is subject sequence name of
the tabular blast and in this case it is scaffolds in a MAG.

Each MAG was renamed with a uniq MAG ID append to beginning of each
scaffold sequence name.

ex: metabat_EN_58_trim_36_scaffold_246

This script splits by "_" and grabs 1-5.

ex: metabat_EN_58_trim_36

New files would {prefix}-metabat_EN_58_trim_36.fltrdBstHts.blst

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

def de_concatenate_tabblast(infile, sampleID):

    # dict of unique column ids from file 1 to look for in file 2
    data = defaultdict(list)

    with open(infile, 'r') as file:
        for line in file:
            X = line.split('\t')
            name = X[1].split('_')
            uid = '_'.join(name[0:5])
            data[uid].append(line)

    for uid, results in data.items():
        Path(uid).mkdir(exist_ok=True)
        outfilename = f'{uid}/{sampleID}-{uid}.fltrdBstHts.blast'
        with open(outfilename, 'w') as outfile:
            for result in results:
                outfile.write(result)

def main():

        # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the tabular blast file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--metagenome_sample_ID',
        help='Please specify the metagenome sample ID!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    de_concatenate_tabblast(
                    args['input_file'],
                    args['metagenome_sample_ID'],
                    )

if __name__ == "__main__":
    main()