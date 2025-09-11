#!/usr/bin/env python

''' Split coupled paired read file in pair 1 and pair 2 files

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: March 31st, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse

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

def Fasta_De_interlace_reads(infile, prefix):

    count = 1
    P1name = prefix + '_P1.fasta'
    P2name = prefix + '_P2.fasta'

    with open(infile, 'r') as f1:
        with open(P1name, 'w') as P1out:
            with open(P2name, 'w') as P2out:
                for name, seq in read_fasta(f1):
                    if count % 2 == 0:
                        P2out.write(name + '\n' + seq + '\n')
                        count += 1
                    else:
                        P1out.write(name + '\n' + seq + '\n')
                        count += 1

def main():
    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the coupled paired read fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the prefix to use for the output files!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    infile = args['input_file']
    prefix = args['output_prefix']
    Fasta_De_interlace_reads(infile, prefix)
    
if __name__ == "__main__":
    main()

