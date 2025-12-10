#!/usr/bin/env python

''' Quick Match Columns of 2 tsv files.

For two tsv files:
Select lines from file 2 that match a column in file 1.
Print matching lines from file 2 to output file.

Takes 7 positional arguments.
file1, column number, header
file2, column number, header
output file

* Column numbers start at 0 (ie: python is 0 indexed)

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 21st, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse

def match_columns(f1, c1, h1, f2, c2, h2, outfile):

    # dict of unique column ids from file 1 to look for in file 2
    match = {}

    # Retrieve unique ids from f1, c1 and populate match dict
    with open(f1, 'r') as f:
        if h1 == 1: header = f.readline()
        for line in f:
            match_id = line.rstrip().split('\t')[c1]
            match[match_id] = ''

    # Look for match dict keys in file 2 and write new file
    with open(f2, 'r') as f, open(outfile, 'w') as o:
        if h2 == 1:
            header = f.readline()
            o.write(header)
        for line in f:
            X = line.rstrip().split('\t')[c2]
            if X in match:
                o.write(line)

def main():

        # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-f1', '--input_file_1',
        help='Please specify the file to retrieve unique IDs from!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c1', '--select_column_1',
        help='Please specify the column number to select (0 indexed)!',
        metavar=':',
        type=int,
        required=False,
        default=0
        )
    parser.add_argument(
        '-h1', '--file_1_header',
        help='File 1 header? (default 0 for no. Use 1 for yes)',
        metavar=':',
        type=int,
        required=False,
        default=0
        )
    parser.add_argument(
        '-f2', '--input_file_2',
        help='Please specify the file to match unique IDs to!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c2', '--select_column_2',
        help='Please specify the column number to match (0 indexed)!',
        metavar=':',
        type=int,
        required=False,
        default=0
        )
    parser.add_argument(
        '-h2', '--file_2_header',
        help='File 2 header? (default 0 for no. Use 1 for yes)',
        metavar=':',
        type=int,
        required=False,
        default=0
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    match_columns(
                    args['input_file_1'],
                    args['select_column_1'],
                    args['file_1_header'],
                    args['input_file_2'],
                    args['select_column_2'],
                    args['file_2_header'],
                    args['output_file_name']
                    )

if __name__ == "__main__":
    main()