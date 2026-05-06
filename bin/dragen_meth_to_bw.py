#!/usr/bin/env python3

import pandas as pd
import pyranges as pr
import pyBigWig
import argparse
import gzip
import os, sys

# add version
__version__ = "0.1.0"


def read_chrom_sizes(fai_file):
    chrom_sizes = []
    if not fai_file.endswith("fai"):
        if os.path.exists(fai_file + ".fai"):
            fai_file = fai_file + ".fai"
        else:
            print(f"Error: Chromosome sizes file '{fai_file}' not found.")
            sys.exit(1)

    with open(fai_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[0]
            size = int(fields[1])
            chrom_sizes.append((chrom, size))
    return chrom_sizes


def add_to_bw(bw, df, coverage=False):
    # Convert numeric columns to appropriate types
    df[["Start", "meth", "unmeth"]] = df[["Start", "meth", "unmeth"]].apply(pd.to_numeric)
    # Calculate start and end columns
    df["End"] = df["Start"]
    df["Start"] = df["Start"] - 1

    # Aggregate rows
    grouped = df.groupby(["Chromosome", "Start", "End"]).agg({"meth": "sum", "unmeth": "sum"}).reset_index()
    # Calculate the required ratio
    grouped["ratio"] = grouped["meth"] / (grouped["meth"] + grouped["unmeth"])
    grouped["total"] = grouped["meth"] + grouped["unmeth"]

    gr = pr.PyRanges(grouped)
    sorted_gr = gr.sort()
    values = sorted_gr.df["total"].tolist() if coverage else sorted_gr.df["ratio"].tolist()
    bw.addEntries(
        sorted_gr.df["Chromosome"].tolist(),
        sorted_gr.df["Start"].tolist(),
        ends=sorted_gr.df["End"].tolist(),
        values=values,
    )


def main():
    parser = argparse.ArgumentParser(description="Process a file and output aggregated data.")
    parser.add_argument(
        "-c", "--coverage", action="store_true", help="Generate coverage bigwig in addition to methylation"
    )
    parser.add_argument("-s", "--chromsizes", type=str, required=True, help="FAI file with chromosome sizes")
    parser.add_argument("-o", "--output_file", type=str, default=None, help="Output file path (default: stdout)")
    # add version
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("input_file", help="Input file path")

    args = parser.parse_args()

    # Read chromosome sizes from FAI file
    chrom_sizes = read_chrom_sizes(args.chromsizes)

    # Read the file into a DataFrame
    cols = ["Chromosome", "Start", "Strand", "meth", "unmeth", "type", "context"]

    # make bigwig file
    bw = pyBigWig.open(args.output_file, "w")
    # Add header information from the chrom sizes file
    bw.addHeader(chrom_sizes)

    covbw = None
    if args.coverage:
        covbw = pyBigWig.open(args.output_file.replace(".bw", ".coverage.bw"), "w")
        covbw.addHeader(chrom_sizes)

    # Initialize the list of rows
    rows = []
    open_func = gzip.open if args.input_file.endswith(".gz") else open

    # Read the file line by line
    with open_func(args.input_file, "rt") as file:
        for line in file:

            data = line.strip().split("\t")

            if len(rows) > 0 and data[0] != rows[-1][0]:  # if new chromosome
                # print chrom to stderr
                print(f"Writing chromosome {rows[-1][0]}", file=sys.stderr)
                add_to_bw(bw, pd.DataFrame(rows, columns=cols))
                if covbw:
                    add_to_bw(covbw, pd.DataFrame(rows, columns=cols), coverage=True)

                rows = []

            if data[5] == "CG":  # Check if 'type' is 'CG'
                # Append the relevant data to the DataFrame
                rows.append(data)

    # Add the last chromosome
    if len(rows) > 0:
        add_to_bw(bw, pd.DataFrame(rows, columns=cols))
        if covbw:
            add_to_bw(covbw, pd.DataFrame(rows, columns=cols), coverage=True)

    bw.close()
    if covbw:
        covbw.close()


if __name__ == "__main__":
    main()
