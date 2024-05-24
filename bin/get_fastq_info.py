#!/usr/bin/env python3

import sys
import os
import csv
import gzip

def extract_info_from_fastq(fastq_file):
    with gzip.open(fastq_file, 'rt') as f:
        first_line = f.readline().strip()
        parts = first_line.split(':')
        flowcell_id = parts[2]
        lane_number = parts[3]
        index_sequence = parts[-1].split('+')
        index1 = index_sequence[0]
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        seq = index_sequence[1]
        reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
        index = index1 + '-' + reverse_complement

    return flowcell_id, lane_number, index

def main():
    sample_id = sys.argv[1]
    fastq1 = sys.argv[2]
    fastq2 = sys.argv[3]

    flowcell_id_1, lane_number_1, index_sequence_1 = extract_info_from_fastq(fastq1)
    flowcell_id_2, lane_number_2, index_sequence_2 = extract_info_from_fastq(fastq2)

    if flowcell_id_1 != flowcell_id_2 or lane_number_1 != lane_number_2:
        print("Error: FASTQ files have different flowcell ID or lane number.")
        sys.exit(1)

    output_file = "fastq_info.csv"
    header = ["sample_id","flowcell_id", "lane_number", "index_sequence", "read1", "read2"]
    data = [[sample_id, flowcell_id_1, lane_number_1, index_sequence_1, os.path.abspath(fastq1), os.path.abspath(fastq2)]]

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerows(data)

    print(f"FASTQ information written to {output_file}")

if __name__ == "__main__":
    main()
