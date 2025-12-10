#!/usr/local/pacerepov1/python/2.7/bin/python

## USAGE :: python scriptname.py file.fasta prefix
## Reads fasta file and replaces sequence names with prefix_# counting from 1 to the end.
## Writes file.ref matching original names to new names.

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

def Filter_prodigal_fasta(infile, minlen):

    X = infile.split('.')

    passedout = X[0] + '_passed.' + X[-1]
    failedout = X[0] + '_failed.' + X[-1]
    po = open(passedout, 'w')
    fo = open(failedout, 'w')

    count_total = 0
    count_pass = 0

    with open(infile, 'r') as f:

        for name, seq in read_fasta(f):
            count_total += 1
            partial = name.split('partial=')[1].split(';')[0]
            length = len(seq)

            lineout = f'{name}\n{seq}\n'

            if length >= minlen and partial != '11':
                count_pass += 1
                po.write(lineout)

            else: fo.write(lineout)

    print(
        f'\n\n\t\tGenes Total: {count_total}\n'
        f'\t\tGenes Passed: {count_pass}\n'
        f'\t\tGenes Filtered: {count_total - count_pass}\n\n\n'
        )

    po.close()
    fo.close()


def main():
    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--fasta_input_file',
        help='Please specify the fasta input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--minimum_sequence_length',
        help='Please specify the minimum sequence length (default=300)!',
        metavar=':',
        type=int,
        default=300,
        required=False
        )
    args=vars(parser.parse_args())

    infile = args['fasta_input_file']
    minlen = args['minimum_sequence_length']

    Filter_prodigal_fasta(infile, minlen)
    
if __name__ == "__main__":
    main()

