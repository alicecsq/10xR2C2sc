#!/usr/bin/env python3

import os
import sys
import argparse

def argParser():
    '''Parses arguments'''
    parser = argparse.ArgumentParser(description = 'Copies the cellranger output for use with Seurat.',
                                     add_help = True,
                                     prefix_chars = '-')
    #parser.add_argument('-a', type=str, action='store',
            #help='Gencode annotation file (gtf). eg. gencode.v29.annotation.gtf')
    #parser.add_argument('-b', type=str, action='store',
            #help='bcGuide file produced from demuxing.')
    parser.add_argument('-e', type=str, action='store',
            help='Output file from featureCounts.')
    #parser.add_argument('-o', type=str, action='store',
            #help='Output path. eg. hg38')
    return vars(parser.parse_args())

args = argParser()
#anno = args['a'] #gtf annotation file
#bcGuide = args['b'] #bcGuide file of all the x most frequent barcodes found from the detBarcodes.py script
exp = args['e'] #featurecounts output file
#outPath = args['o']
#if not outPath:
    #outPath = ""
#if not outPath.endswith('/'):
    #outPath += '/'

def modifyFC(fcIn): #this basically gets the counts per gene (we want a way to store the cell number info for each count column too I would think)
    #We can make a second dictionary that holds all the cell names and the index of cell name in this dict vs the gene count id ones can be the same. 
    '''
    Reads the featureCounts output and throws it into a dictionary
    countDict = {geneID: [0, 1, 0, 0, ...], ...}
    The index of the list +1 is the cell # #we can just have the cell number be the column header because that is what it is in mine
    '''
    first = True
    countDict = {}
    sys.stderr.write('Reading in the featureCounts output\n')
    for line in open(fcIn): #for each line in the featurecounts file...
        if line.startswith('#'): # featureCounts command info
            continue
        if first:
            first = False
            cellnames=line.rstrip().split('\t')[6:]
            cellnumindex=[]
            for i in range(0,len(line.rstrip().split('\t')[6:]),1):
                cellnumber=((cellnames[i].split('/')[1]).split('.')[0]).split('_')[1]
                cellnumindex.append(cellnumber) #creates cellnum list that has the cell number in the same index position as the genecount dict entries
            continue
        print(cellnumindex)
        line = line.rstrip().split('\t')
        gene = line[0]
        counts = line[6:] # only need the count columns
        countDict[gene] = counts #creates dict with gene:count/cell key:value pairs. Entries have same indexing as cellnumindex so counts can be mapped back to the right cell
    return countDict, cellnumindex


countDict=modifyFC(exp)
#for entry in countDict:
    #print(len(countDict[entry]))
