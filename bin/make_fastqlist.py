#!/usr/bin/env python3

import xml.etree.ElementTree as ET
import json, sys, os, csv, argparse, gzip, re
from collections import Counter
import pandas as pd

__version__ = '1.0.0'

"""
This script makes a custom-formatted fastq_list.csv file from either an existing fastq_list.csv,
a dragen demux directory, or read1 and read2 FASTQ files. Features include:

-Accepts a fastq_list.csv file in this format:

RGID,RGSM,RGLB,Lane,Read1File,Read2File

-Returns a fastq_list.csv file in this format:

RGID,RGSM,RGLB,Lane,RGPL,Read1File,Read2File

-Generates tags with this format:

RGID: <flowcell>.<i7index>.<i5index>.<lane>
RGSM: <sample id> (passed as argument)
RGLB: <sample id>.<i7index>.<i5index>
Lane: <lane number>
RGL: <runid>.<instrument><side>.<flowcell>.<flowcelltype>.<flowcell_lot>.<reagent_lot>.<runrecipe>
Read1File: <read1 file path>
Read2File: <read2 file path>

-Obtains RGPL info from either a RunParameters.xml file or the reads (info from reads is incomplete)
-If demux path is given, look for Reports/fastq_list.csv

"""


def checkfile(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    return file_path

# NOTE(dhs): get flowcell, lane, read length from reads to make runinfo
def make_runinfo_from_read(readpath):

    # open read 1 fastq and read header and get seq length
    readName = ''
    with gzip.open(readpath, 'rt') as gz_file:
        readName = gz_file.readline().strip()
    
    # if indexes are present in the read header, get the first 100
    # and find the most common one (to account for mismatches/errors in index read)
    parts = readName.split(':')
    index1 = 'UNKNOWN'
    index2 = 'UNKNOWN'
    index1len = '?'
    index2len = '?'
    indexes = readName.split(' ')
    readlen = '?'
    indexlist = []
    seqlist = []
    with gzip.open(readpath, 'rt') as file:
        for i, line in enumerate(file):
            if i % 4 == 0:  # Read names are on every 4th line starting from 0
                read_name = line.strip()
                indexes = read_name.split(' ')
                if len(indexes) > 1:
                    index = indexes[1].split(':')[-1]
                    indexlist = indexlist + [index]
    
            if i > 0 and i % 1 == 0:
                seqlist.append(len(line.strip()))

            if i == 399:  # Stop after reading the first 100 entries (0 to 399)
                break

    if len(indexlist) > 0:
        counter = Counter(indexlist)
        most_common_index, _ = counter.most_common(1)[0]
        index1, index2 = most_common_index.split('+')

    counter = Counter(seqlist)
    readlen, _ = counter.most_common(1)[0]

    # if indexes are ?, try to recover from the read name
    if index1 == 'UNKNOWN':
        pattern = '[ACTG]{8,}'
        matches = re.findall(pattern, os.path.basename(readpath))
        if len(matches) == 2:
            index1, index2 = matches
            index1len = f"{len(index1)}?"
            index2len = f"{len(index2)}?"
        elif len(matches) == 1:
            index1 = matches
            index2 = matches
        else:
            index1 = f"{os.path.basename(readpath).split('.')[0]}.{index1}"

    # Make run info dict
    runinfo = {'RunId':f"RUN_{parts[0][1:]}_{str(int(parts[1])).zfill(4)}_{parts[2]}",
                'Flowcell':parts[2],
                'Lane':parts[3],
                'Instrument':parts[0][1:],
                'Read1Cycles':f"{readlen}",
                "Index1Cycles": index1len,
                "Index1Reverse": "?",
                "Index2Cycles": index2len,
                "Index2Reverse": "?",
                "Read2Cycles": f"{readlen}",
                "FlowcellType":  "?",
                "FlowcellLot": "?",
                "ReagentLot": "?",
                "Index1": index1,
                "Index2": index2
    }

    return runinfo

# NOTE(dhs): Make run info dict from files in demux path
def make_runinfo_from_run_parameters(filepath):

    runinfo = {}
    
    tree = ET.parse(filepath)
    root = tree.getroot()
    runinfo = {}
    runinfo['RunId'] = root.find('RunId').text
    for el in root.findall('.//PlannedReads')[0].findall('Read'):
        runinfo[el.attrib['ReadName'] + 'Cycles'] = int(el.attrib['Cycles'])

    runinfo['FlowcellType'] = root.find('RecipeName').text.split(' ')[0]
    runinfo['Instrument'] = root.find('InstrumentSerialNumber').text 
    runinfo['Side'] = root.find('Side').text
    # fine the serial number of the flowcell:
    for el in root.findall('ConsumableInfo')[0].findall('ConsumableInfo'):
        el2 = el.find('Type')
        if el2 is not None:
            if el2.text == 'FlowCell':
                runinfo['Flowcell'] = el.find('SerialNumber').text
                runinfo['FlowcellLot'] = el.find('LotNumber').text
            elif el2.text == 'Reagent':
                runinfo['ReagentLot'] = el.find('LotNumber').text
    
    runinfo['Instrument'] = ''.join([runinfo.get('Instrument') or '?',runinfo.get('Side') or '?'])

    return runinfo



def main():

    parser = argparse.ArgumentParser(description="Prepare fastq_list file from a passed list file, Dragen demux directory, or reads")
    parser.add_argument('-i', '--id', type=str, required=True, help="Sample ID")
    parser.add_argument('-1','--read1', type=checkfile, help="Path to read1")
    parser.add_argument('-2','--read2', type=checkfile, help="Path to read2")
    parser.add_argument('-f','--fastqlist', type=checkfile, help="Path to fastq_list.csv")
    parser.add_argument('-d','--demuxpath', type=str, help="Path to Dragen demux directory")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)

    args = parser.parse_args()

    id = args.id

    # option 1, if read1 and read2 are supplied
    if args.read1 and args.read2:

        read1 = os.path.realpath(args.read1)
        read2 = os.path.realpath(args.read2)

        runinfo = make_runinfo_from_read(read1)
        runinfo2 = make_runinfo_from_read(read2)
        
        if (runinfo['Lane'] != runinfo2['Lane'] or 
            runinfo['Flowcell'] != runinfo2['Flowcell'] or
            runinfo['Index1'] != runinfo2['Index1']):
            print(f"Read 1 and 2 do not match\n\t{read1}\n\t{read2}")
            sys.exit(1)

        fqlistout = pd.DataFrame(columns=['RGID','RGSM','RGLB','Lane','RGPL','Read1File','Read2File'])

        rgid = '.'.join([runinfo['Flowcell'],runinfo['Index1'],runinfo['Index2'],runinfo['Lane']])
        rglb = '.'.join([id,runinfo['Index1'],runinfo['Index2']])
        rgpl = f"{runinfo['RunId']}.{runinfo['Instrument']}.{runinfo['Flowcell']}.{runinfo['FlowcellType']}.{runinfo['FlowcellLot']}.{runinfo['ReagentLot']}.{runinfo['Read1Cycles']}x{runinfo['Index1Cycles']}x{runinfo['Index2Cycles']}x{runinfo['Read2Cycles']}"

        fqlistout.loc[len(fqlistout)] = [rgid,id,rglb,runinfo['Lane'],rgpl,read1,read2]

        fqlistout.to_csv(f'{id}.fastq_list.csv',index=False)

    # options 2 and 3, if a fastq list or demux path are supplied
    else:
        
        runparams_path = None
        fastqlistFile = None
        runinfo = {}
        
        if args.demuxpath:
            # locate the RunParameters.xml file and check if it exists:
            runparams_path = os.path.join(args.demuxpath, 'Reports','RunParameters.xml')
            if not os.path.exists(runparams_path):
                raise ValueError('fastq_list.csv file not found in the specified directory')

            fastqlistFile = os.path.join(args.demuxpath,'Reports','fastq_list.csv')

        elif args.fastqlist:
            fastqlistFile = args.fastqlist

        if not os.path.exists(fastqlistFile):
            raise ValueError('fastq_list.csv file not found in the specified directory')

        if runparams_path and os.path.exists(runparams_path):
            runinfo = make_runinfo_from_run_parameters(runparams_path)

        fqlist = pd.read_csv(fastqlistFile)

        fqlist = fqlist[fqlist['RGSM']==id]

        fqlistout = pd.DataFrame(columns=['RGID','RGSM','RGLB','Lane','RGPL','Read1File','Read2File'])

        for index, row in fqlist.iterrows():
            read1 = row['Read1File']
            read2 = row['Read2File']
            lane = row['Lane']

            if not os.path.exists(read1):
                if os.path.exists(os.path.join(args.demuxpath,os.path.basename(read1))):
                    read1 = os.path.join(args.demuxpath,os.path.basename(read1))
                else:
                    print(f"Read 1 file: {os.path.basename(read1)} does not exist")
                    sys.exit(1)
            
            if not os.path.exists(read2):
                if os.path.exists(os.path.join(args.demuxpath,os.path.basename(read2))):
                    read2 = os.path.join(args.demuxpath,os.path.basename(read2))
                else:
                    print(f"Read 2 file: {os.path.basename(read2)} does not exist")
                    sys.exit(1) 

            if not runinfo:
                runinfo = make_runinfo_from_read(read1)

            rgid = f"{runinfo['Flowcell']}.{row['RGID']}"
            rglb = f"{id}.{row['RGID'][:-2]}"
            rgpl = f"{runinfo['RunId']}.{runinfo['Instrument']}.{runinfo['Flowcell']}.{runinfo['FlowcellType']}.{runinfo['FlowcellLot']}.{runinfo['ReagentLot']}.{runinfo['Read1Cycles']}x{runinfo['Index1Cycles']}x{runinfo['Index2Cycles']}x{runinfo['Read2Cycles']}"

            fqlistout.loc[len(fqlistout)] = [rgid,id,rglb,lane,rgpl,read1,read2]

        fqlistout.to_csv(f'{id}.fastq_list.csv',index=False)

if __name__ == "__main__":
    main()

