#!/usr/bin/env python3

import argparse
import os
import sys

import pandas as pd

__version__ = "1.0.1"

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process spreadsheet file for samples that require demultiplexing, alignment, or summary analysis."
    )
    parser.add_argument(
        "--input_file", "-i", required=True, help="Path to the input CSV file"
    )
    parser.add_argument(
        "--output_dir", "-o", help="Directory to save the output CSV files"
    )
    return parser.parse_args()


def clean_sample_name(value: str) -> str:
    """Remove whitespace and replace spaces with underscores."""
    if isinstance(value, str):
        return value.strip().replace(" ", "_")
    return value


def save_df(dataframe: pd.DataFrame, output_dir: str, filename: str) -> None:
    """Save dataframe to a CSV file, applying a cleaning function if provided."""
    if not dataframe.empty:
        dataframe.loc[:, "id"] = dataframe["id"].apply(
            clean_sample_name
        )

        dataframe = dataframe.dropna(axis=1, how="all")
        dataframe = dataframe.loc[:, dataframe.ne("").any()]

        dataframe.to_csv(f"{output_dir}/{filename}", index=False)


def get_read_pairs_from_samplemap(samplemap_path: str, match_id: str) -> list:
    """Return list of (read1, read2) absolute path tuples from a samplemap for a given Library Name."""
    smap = pd.read_csv(samplemap_path)
    smap_dir = os.path.realpath(os.path.dirname(samplemap_path))
    matched = smap[smap['Library Name'] == match_id]
    pairs = []

    if 'FASTQ Path - Read 1' in smap.columns and 'FASTQ Path - Read 2' in smap.columns:
        for _, row in matched.iterrows():
            pairs.append((
                os.path.join(smap_dir, row['FASTQ Path - Read 1']),
                os.path.join(smap_dir, row['FASTQ Path - Read 2']),
            ))
    elif 'FASTQ' in smap.columns:
        r1_rows = matched[matched['FASTQ'].str.contains('_R1_', na=False)]
        r2_rows = matched[matched['FASTQ'].str.contains('_R2_', na=False)]
        key_cols = [c for c in ['Flowcell ID', 'Index Sequence', 'Flowcell Lane'] if c in smap.columns]
        for _, r1_row in r1_rows.iterrows():
            r2_match = r2_rows
            for col in key_cols:
                r2_match = r2_match[r2_match[col] == r1_row[col]]
            if not r2_match.empty:
                pairs.append((
                    os.path.join(smap_dir, r1_row['FASTQ']),
                    os.path.join(smap_dir, r2_match.iloc[0]['FASTQ']),
                ))

    return pairs


def process_input(input_file: str, output_dir: str) -> None:
    """Process input file, merge with sample info, and create output files."""
    ext = input_file.split(".")[-1]
    df = pd.read_csv(input_file, sep="\t" if ext == "tsv" else ",",engine="python",on_bad_lines="error")
            
    # modify sex column
    if "sex" in df.columns:
        df["sex"] = df["sex"].str.lower().replace({"m": "male", "f": "female"})
        df.dropna(axis=1, how="all", inplace=True)

    # Pivot from wide tumor/normal format to long format
    if "tumor_id" in df.columns and "normal_id" in df.columns:
        common_cols = [c for c in df.columns if not c.startswith('tumor_') and not c.startswith('normal_')]
        
        tumor_df = df[common_cols].copy()
        tumor_df['sample_type'] = 'tumor'
        normal_df = df[common_cols].copy()
        normal_df['sample_type'] = 'normal'
        
        tumor_df['id'] = df['tumor_id']
        normal_df['id'] = df['normal_id']
        for base_name in ['read1', 'read2', 'fastq_list', 'cram', 'bam']:
            if f'tumor_{base_name}' in df.columns and f'normal_{base_name}' in df.columns:
                tumor_df[base_name] = df[f'tumor_{base_name}']
                normal_df[base_name] = df[f'normal_{base_name}']

        if 'fastq_list' in df.columns:
            tumor_df['fastq_list'] = df['fastq_list']
            normal_df['fastq_list'] = df['fastq_list']

        elif 'samplemap' in df.columns:
            tumor_df['samplemap'] = df['samplemap']
            normal_df['samplemap'] = df['samplemap']

        df = pd.concat([tumor_df, normal_df]).reset_index(drop=True)

        if 'sample_id' not in df.columns:
            df['sample_id'] = df['id']
            
    # remove columns if not valid header
    valid_headers = [
            "id",
            "individual_id",
            "sample_type",
            "sample_id"
            "sex",
            "read1",
            "read2",
            "fastq_list",
            "samplemap",
            "cram",
            "bam",
            "dragen_path"
    ]

    df = df.loc[:, [col for col in df.columns if col in valid_headers]]

    if "dragen_path" in df.columns:
        analysis_cols = ["id", "dragen_path"]
        analysis_df = df.dropna(
            subset=[col for col in analysis_cols if col in df], thresh=1
        )
        save_df(analysis_df, output_dir, "analysis_samples.csv")

    else:
        # Alignment samples. Note this includes samples to demux as well as already demuxed fastqs and cram/bam for realignment.
        alignment_cols = ["bam", "cram", "read1", "read2", "fastq_list", "samplemap"]
        alignment_df = df.dropna(
            subset=[col for col in alignment_cols if col in df], thresh=1
        )

        if 'samplemap' in alignment_df.columns and alignment_df['samplemap'].notna().any():
            expanded_rows = []
            for _, row in alignment_df.iterrows():
                if pd.isna(row.get('samplemap')):
                    expanded_rows.append(row)
                    continue
                use_sample_id = (
                    'sample_type' in alignment_df.columns
                    and pd.notna(row.get('sample_type'))
                    and 'sample_id' in alignment_df.columns
                    and pd.notna(row.get('sample_id'))
                )
                match_id = row['sample_id'] if use_sample_id else row['id']
                pairs = get_read_pairs_from_samplemap(row['samplemap'], match_id)
                if pairs:
                    for read1, read2 in pairs:
                        new_row = row.copy()
                        new_row['read1'] = read1
                        new_row['read2'] = read2
                        expanded_rows.append(new_row)
                else:
                    sys.exit(f"ERROR: ID '{match_id}' not found in samplemap '{row['samplemap']}'")
            alignment_df = pd.DataFrame(expanded_rows).reset_index(drop=True)
            alignment_df.drop(columns=['samplemap'], inplace=True)

        save_df(alignment_df, output_dir, "alignment_samples.csv")

def main() -> None:
    args = parse_arguments()

    output_dir = args.output_dir
    if not output_dir or not os.path.exists(output_dir):
        output_dir = os.getcwd()

    process_input(args.input_file, output_dir)


if __name__ == "__main__":
    main()
