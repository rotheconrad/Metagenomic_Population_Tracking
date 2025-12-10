#!/usr/bin/env python

'''Calculate Density from Temp and Salinity metadata

Written for GoM Project Metagenome_Samples_MetaData.tsv file.
TARA, HOT, Geotrace sample metadata provide temp and salinity
but not density. Needed a density column for some plotting and analysis.

Writes out new metadata file with density column.
Metagenome_Samples_MetaData_Desnity.tsv

Need to install python seawater or gsw package to calculate density:
https://pypi.org/project/seawater/
conda install -c conda-forge seawater or pip install seawater

https://pypi.org/project/gsw/
conda install -c conda-forge gsw or pip install gsw

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 02 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seawater

def calculate_density(metadata):

    # import metadata file as pandas datafram
    df = pd.read_csv(metadata, sep='\t')

    # Extract the Temp and Salinity columns as numpy arrays
    temp = df.Temp.to_numpy()
    salt = df.Salinity.to_numpy()

    # calculate density from salinity and temp
    rho = seawater.eos80.dens0(salt, temp) # (salinity, temp)

    # Substract 1000 to convert to sigma-t
    rho = rho - 1000

    # append density column to dataframe
    df['Density'] = rho

    return df

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-m', '--metadata_input_file',
        help='Please specify the metadata input file name!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # define input params
    metadata = args['metadata_input_file']

    # Do the density!
    mdf = calculate_density(metadata)

    outname = metadata.split('.')[0] + "_Density.tsv"
    mdf.to_csv(outname, sep='\t', header=True, index=False)

if __name__ == "__main__":
    main()
