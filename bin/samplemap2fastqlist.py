#!/usr/bin/env python3

import argparse
import csv
import os
import sys


def detect_format(fieldnames):
    if "FASTQ Path - Read 1" in fieldnames and "FASTQ Path - Read 2" in fieldnames:
        return "samplemap2"
    elif "FASTQ" in fieldnames:
        return "samplemap"
    else:
        raise ValueError(f"Unrecognized samplemap format. Headers: {fieldnames}")


def process_samplemap2(csv_reader, csv_writer, input_dir, prefix=None):
    for row in csv_reader:
        RGID = f"{row['Flowcell ID']}.{row['Index Sequence']}.{row['Flowcell Lane']}"
        RGSM = row["Library Name"]
        if prefix:
            RGSM = prefix + "." + RGSM
        RGLB = f"{row['Library Name']}.{row['Index Sequence']}"
        Lane = row["Flowcell Lane"]
        Read1File = os.path.join(input_dir, row["FASTQ Path - Read 1"])
        Read2File = os.path.join(input_dir, row["FASTQ Path - Read 2"])
        csv_writer.writerow(
            {"RGID": RGID, "RGSM": RGSM, "RGLB": RGLB, "Lane": Lane, "Read1File": Read1File, "Read2File": Read2File}
        )


def process_samplemap(csv_reader, csv_writer, input_dir, prefix=None):
    groups = {}
    for row in csv_reader:
        key = (row["Flowcell ID"], row["Index Sequence"], row["Flowcell Lane"], row["Library Name"])
        fastq = row["FASTQ"]
        if "_R1_" in fastq:
            groups.setdefault(key, {})["R1"] = fastq
        elif "_R2_" in fastq:
            groups.setdefault(key, {})["R2"] = fastq

    for key, files in groups.items():
        flowcell_id, index_seq, lane, lib_name = key
        if "R1" not in files or "R2" not in files:
            print(f"Warning: missing R1 or R2 for {key}", file=sys.stderr)
            continue
        RGID = f"{flowcell_id}.{index_seq}.{lane}"
        RGSM = lib_name
        if prefix:
            RGSM = prefix + "." + RGSM
        RGLB = f"{lib_name}.{index_seq}"
        Read1File = os.path.join(input_dir, files["R1"])
        Read2File = os.path.join(input_dir, files["R2"])
        csv_writer.writerow(
            {"RGID": RGID, "RGSM": RGSM, "RGLB": RGLB, "Lane": lane, "Read1File": Read1File, "Read2File": Read2File}
        )


def process_input_file(input_file, prefix=None):
    with open(input_file, mode="r") as file:
        csv_reader = csv.DictReader(file)
        fmt = detect_format(csv_reader.fieldnames)

        fieldnames = ["RGID", "RGSM", "RGLB", "Lane", "Read1File", "Read2File"]
        csv_writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
        csv_writer.writeheader()

        input_dir = os.path.realpath(os.path.dirname(input_file))

        if fmt == "samplemap2":
            process_samplemap2(csv_reader, csv_writer, input_dir, prefix)
        else:
            process_samplemap(csv_reader, csv_writer, input_dir, prefix)


def main():
    parser = argparse.ArgumentParser(description="Process a FASTQ metadata file and output CSV.")
    parser.add_argument("input_file", type=str, help="Path to the input file")
    parser.add_argument("-p", "--prefix", type=str, help="Prefix for library name.")
    args = parser.parse_args()

    process_input_file(args.input_file, args.prefix)


if __name__ == "__main__":
    main()
