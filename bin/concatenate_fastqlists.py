#!/usr/bin/env python

import os, sys, csv

if len(sys.argv) < 3:
    print("Usage: concatenate_fastqlists.py <meta.id> <file1.csv> <file2.csv> ...")
    sys.exit(1)

id = sys.argv[1]
csv_files = sys.argv[2:]

# Define output file names
output_csv = f"{id}.fastq_list.csv"
output_read_files = f"{id}.fastqs.csv"

# Prepare the header for the first output file
header = ["RGID", "RGSM", "RGLB", "Lane", "RGPL", "Read1File", "Read2File"]

output_rows = []
read_files = []

# Process each CSV file
for csv_file in csv_files:
    with open(csv_file, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 7:
                read_files = read_files + [row[5], row[6]]
                row[5] = os.path.basename(row[5])
                row[6] = os.path.basename(row[6])
                output_rows.append(row)

# Write to the first output file with the specified header
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(output_rows)

# Write to the second output file with Read1File and Read2File, newline delimited
with open(output_read_files, "w") as f:
    for read_file in read_files:
        f.write(read_file + "\n")

