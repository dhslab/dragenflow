#!/usr/bin/env python3

# write a script using argparse and pandas that accepts a table of transcripts with this format: 
# gene_id gene_name       transcript_id   transcript_name
# and a table of gene expression values like this: 
# Name    Length  EffectiveLength TPM     NumReads
# and then joins the two tables on the gene_id column and writes the result to a new file.

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", type=str, help="expression file")
    parser.add_argument("-d","--database", type=str, help="transcript database flatfile")
    parser.add_argument("-t","--transcripts", action='store_true',default=False, help="Convert transcripts instead of genes")
    parser.add_argument("-o","--outfile", help="output file")
    args = parser.parse_args()
    
    transcripts = pd.DataFrame()
    if args.transcripts:
        transcripts = pd.read_csv(args.database,sep="\t")[["transcript_id","gene_id","gene_name","transcript_name"]].drop_duplicates()
        trx_expression = pd.read_csv(args.input,sep="\t")
        outdf = trx_expression.merge(transcripts,left_on="Name",right_on="transcript_id",how="left")
        outdf.loc[outdf['transcript_name'].isna(),"transcript_name"] = outdf.loc[outdf['transcript_name'].isna(),"transcript_id"]
        outdf.rename(columns={'Name':'Id','transcript_name':'Name','gene_id':'GeneId','gene_name':'GeneName'},inplace=True)
        outdf["Id Name GeneId GeneName Length EffectiveLength TPM NumReads".split(" ")].to_csv(args.outfile,sep="\t",index=False)

    else:
        transcripts = pd.read_csv(args.database,sep="\t")[["gene_id","gene_name"]].drop_duplicates()
        gene_expression = pd.read_csv(args.input,sep="\t")
        outdf = gene_expression.merge(transcripts,left_on="Name",right_on="gene_id",how="left")
        outdf.loc[outdf['gene_name'].isna(),"gene_name"] = outdf.loc[outdf['gene_name'].isna(),"gene_id"]
        outdf.rename(columns={'Name':'Id','gene_name':'Name'},inplace=True)
        outdf["Id Name Length EffectiveLength TPM NumReads".split(" ")].to_csv(args.outfile,sep="\t",index=False)

if __name__ == "__main__":
    main()

