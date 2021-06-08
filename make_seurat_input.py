#!/usr/bin/env python3

'''creates barcodes.tsv, features.tsv and matrix file for Seurat input'''

import os
import sys
import argparse

def argParser():
    '''Parses arguments'''
    parser = argparse.ArgumentParser(description = 'Copies the cellranger output for use with Seurat.',
                                     add_help = True,
                                     prefix_chars = '-')
    parser.add_argument('-a', type=str, action='store',
            help='Gencode annotation file (gtf). eg. gencode.v29.annotation.gtf')
    parser.add_argument('-b', type=str, action='store',
            help='bcGuide file produced from demuxing.')
    parser.add_argument('-e', type=str, action='store',
            help='Output file from featureCounts.')
    parser.add_argument('-o', type=str, action='store',
            help='Output path. eg. hg38')
    return vars(parser.parse_args())

args = argParser()
anno = args['a'] #gtf annotation file
bcGuide = args['b'] #bcGuide file of all the x most frequent barcodes found from the detBarcodes.py script
exp = args['e'] #featurecounts output file
outPath = args['o']
if not outPath:
    outPath = ""
if not outPath.endswith('/'):
    outPath += '/'


def makeGenes(inFile): #makes dictionary with genename:linenum key:value pairs
    '''
    Makes gene dictionary out of annotation file and creates 'features.tsv' file with geneID and name. Gene dictionary stores the line number
    in features.tsv where each gene ID is for matrix file formatting.
    '''
    idDict, geneDict = {}, {}
    for line in open(inFile):
        if line[0] == '#':
            continue
        line = line.rstrip().split('\t')
        line = [x.strip() for x in line[8].split(';')]
        id, name = '', ''
        for section in line:
            if section[:7] == 'gene_id':
                id = section.split('"')[1]
            if section[:9] == 'gene_name':
                name = section.split('"')[1]
            if id and name:
                break
        idDict[id] = name

    genesOut = open(outPath + 'features.tsv', 'w+')
    index = 1
    sys.stderr.write('Writing genes to {}features.tsv\n'.format(outPath))
    for gID, gName in idDict.items():
        genesOut.write(gID + '\t' + gName+ '\t'+'Gene Expression'+ '\n') #reflects the new v3.1 features.tsv format
        geneDict[gID] = index
        index += 1
    genesOut.close()
    return geneDict

def makeBC(inFile):
    '''
    Makes the barcodes.tsv file by adding '-1' to the end of each barcode--these numbers change if multiple 10x wells were used!
    '''
    numCells = 0
    bcOut = open(outPath + 'barcodes.tsv', 'w+')
    sys.stderr.write('Writing barcodes to {}barcodes.tsv\n'.format(outPath))
    for line in open(inFile):
        line = line.rstrip().split('\t')
        bcOut.write(line[1] + '-1\n')
        numCells += 1
    bcOut.close()
    return numCells


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
        #print(cellnumindex)
        line = line.rstrip().split('\t')
        gene = line[0]
        counts = line[6:] # only need the count columns
        countDict[gene] = counts #creates dict with gene:count/cell key:value pairs. Entries have same indexing as cellnumindex so counts can be mapped back to the right cell
    return countDict, cellnumindex

def make_mtx(cellnumindex,countDict,geneDict):
    mtxOut = open(outPath + 'matrix.mtx', 'w+')

    genecounts=0 #number of lines in mtx file
    for i in range(0,len(cellnumindex),1):
        for entry in countDict:
            if countDict[entry][i] != '0':
                genecounts += 1
                #print(geneDict[entry],i+1,countDict[entry][i])
    print(len(geneDict)) #number of cells
    mtxOut = open(outPath + 'matrix.mtx', 'w+')
    sys.stderr.write('Writing matrix to {}matrix.mtx\n'.format(outPath))
    mtxOut.write("%%MatrixMarket matrix coordinate integer general\n%\n")
    mtxOut.write("{0} {1} {2}\n".format(str(len(geneDict)), str(len(cellnumindex)), str(genecounts)))

    for i in range(0,len(cellnumindex),1):
        for entry in countDict:
            if countDict[entry][i] != '0':
                mtxOut.write(str(geneDict[entry]) + ' ' + str(i+1) + ' ' + str(countDict[entry][i]) + '\n')
 




makeBC(bcGuide) #make barcodes.tsv file


countDict,cellnumindex=modifyFC(exp) #takes counts per gene per cell info from featurecounts file and puts into dictionary
geneDict = makeGenes(anno) #makes gene dictionary from annotation file and creates features.tsv

make_mtx(cellnumindex,countDict,geneDict) #uses cell number index, counts per gene, and gene dictionary to make MEX file

