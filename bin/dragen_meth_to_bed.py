#!/usr/bin/env python3

import pandas as pd
import pyranges as pr
import argparse
import gzip
import os, sys

# add version 
__version__ = '0.1.0'

def process_file(input_file, output_file):
    # Read the file into a DataFrame
    cols = ['Chromosome', 'pos', 'strand', 'meth', 'unmeth', 'type', 'context']
    df = pd.DataFrame(columns=cols)

    rows = []
    open_func = gzip.open if input_file.endswith('.gz') else open

    # Read the file line by line
    with open_func(input_file, 'rt') as file:
        for line in file:
            data = line.strip().split('\t')
            if data[5] == 'CG':  # Check if 'type' is 'CG'
                # Append the relevant data to the DataFrame
                rows.append(data)

    df = pd.DataFrame(rows, columns=cols)
    # Convert numeric columns to appropriate types
    df[['pos', 'meth', 'unmeth']] = df[['pos', 'meth', 'unmeth']].apply(pd.to_numeric)

    # Calculate start and end columns
    df['Start'] = df['pos'] - 1
    df.loc[df['strand'] == '-', 'Start'] = df['pos'] - 2
    df['End'] = df['Start'] + 2

    # Aggregate rows
    grouped = df.groupby(['Chromosome', 'Start', 'End']).agg({'meth': 'sum', 'unmeth': 'sum'}).reset_index()

    # Calculate the required ratio
    grouped['ratio'] = grouped['meth'] / (grouped['meth'] + grouped['unmeth'])
    grouped['total'] = grouped['meth'] + grouped['unmeth']

    gr = pr.PyRanges(grouped)
    sorted_gr = gr.sort()

    # Determine if output should be gzipped
    if output_file and output_file.endswith('.gz'):
        sorted_gr.df.to_csv(output_file, sep='\t', columns=['Chromosome', 'Start', 'End', 'ratio', 'total'], header=False, index=False, compression='gzip')
    elif output_file:
        sorted_gr.df.to_csv(output_file, sep='\t', columns=['Chromosome', 'Start', 'End', 'ratio', 'total'], header=False, index=False)
    else:
        sorted_gr.df.to_csv(sys.stdout, sep='\t', columns=['Chromosome', 'Start', 'End', 'ratio', 'total'], header=False, index=False)

def main():
    parser = argparse.ArgumentParser(description="Process a file and output aggregated data.")
    parser.add_argument("input_file", help="Input file path")
    parser.add_argument("-o", "--output_file", type=str,default=None, help="Output file path (default: stdout)")
    # add version
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    # Process the file
    process_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

    