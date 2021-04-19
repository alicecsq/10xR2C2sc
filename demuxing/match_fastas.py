#!/usr/bin/env python3
# Roger Volden

'''
Used to match up the fasta entries between the demux fasta and the
full length consensi. This is because when demuxing, sometimes reads
get dropped, which throws off the script that sorts reads into cells.

Assumes the demux fasta file has fewer lines than the flc file.

Usage:
python3 match_fastas.py \
    kmer_demuxed.fasta \
    R2C2_postprocessed.fasta \
    >R2C2_matched.fasta
'''

import sys

def readFasta(inFile): #parses full length consensus read fasta file; basically creates readname:sequence dictionary 
    readDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            lastHead = line[1:].split('_')[0]
            readDict[lastHead] = ''
        else:
            readDict[lastHead] += line
    return readDict

def readDemux(inFile, flcDict):
    for line in open(inFile):
        line = line.rstrip()
        if line.startswith('>'): 
            header = line[1:].split('_')[0] #read rootname
            print('>' + line[1:].split('|')[0]) 
            print(flcDict[header]) #not sure that this is actually the right thing but let's try it...

def main():
    demux, flc = sys.argv[1], sys.argv[2] #demux is the kmer_demuxed.fasta, flc is the full length consensus reads fasta

    flcDict = readFasta(flc)

    readDemux(demux, flcDict)

if __name__ == '__main__':
    main()
