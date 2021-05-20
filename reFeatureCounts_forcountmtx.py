#!/usr/bin/env python3
'''Script for taking UMI-merged reads from Map10xUMIs.py, remapping thesefiltered reads with minimap2 and getting count matrix from FeatureCounts'''
'''dependencies: FeatureCounts, minimap2'''




'''Add FeatureCounts and minimap2 executables to $PATH, otherwise need to put in a config file and unannotate config file fxns'''

import argparse
import editdistance
import numpy as np
import os
import sys
import mappy as mm
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path to individual cell .final.fasta UMI-merged/polished files
parser.add_argument('-o', '--output_path', type=str) #where you want minimap2 and FeatureCounts outputs to go 
parser.add_argument('-g', '--gtf_file', type=str) #same gtf file used for Map10xUMIs.py
parser.add_argument('-c', '--config_file', type=str) #config file containing paths to executables if not in $PATH
parser.add_argument('-r', '--ref_genome', type=str) #reference GENOME for minimap2 alignment (same used for Map10xUMIs.py)


args = parser.parse_args()
input_path=args.input_path
path=args.output_path
config_file=args.config_file
ref_genome=args.ref_genome
gtf=args.gtf_file

'''Parse config'''
def configReader(configIn):
    progs = {}
    for line in open(configIn): 
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t') 
        progs[line[0]] = line[1] 
        
    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'samtools','featureCounts', 'psl2pslx']) #don't think I actually need psl2pslx and emtrey...
    inConfig = set() 
    for key in progs.keys(): 
        inConfig.add(key)
        if key not in possible: 
            raise Exception('Check config file')

    for missing in possible-inConfig: 
        if missing == 'consensus': 
            path = 'consensus.py' 
        else:
            path = missing  
        progs[missing] = path 
        #sys.stderr.write('Using ' + str(missing)
                         #+ ' from your path, not the config file.\n')
    return progs

progs = configReader(config_file)
featureCounts = progs['featureCounts']
minimap2 = progs['minimap2']


'''Parse each cell fasta, trim barcode and UMI sequences from each read and concat into one fasta file'''
def trim (inFile):
    readdict={}
    filename=inFile.split('/')[1]
    cell=filename.split('_')[0]
    number=filename.split('_')[1]
    barcode=filename.split('_')[2]
    bc=barcode.split('.')[0]
    for line in open(inFile):
        if line.startswith('>'):
            name=line.split('_')[0]
            name=name+'_'+cell+'_'+number+'_'+bc
        else:
            sequence=line
            trim=len(sequence)-29
            sequence=sequence[0:trim]
            readdict[name]=sequence
    return readdict


def concat_fasta(input_path, path):
    #sys.stderr.write('Concatenating reads')
    final = open(path + '/allfinalreads.fasta', 'w')
    fileList=[]
    for file in os.listdir(input_path):
        if '.final.fasta' in file:
            fileList.append(file)
    #print(sorted(fileList, key=lambda x: int(x.split('_')[1])))
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])):
        fasta_file = input_path + '/' + file
        readdict=trim(fasta_file)
        for i in readdict:
            sequence=readdict[i] 
            final.write('%s\n%s\n' % (i,sequence))
    return path+'/allfinalreads.fasta'



'''Map reads in the concatenated fasta file with minimap2--uses same parameters as Map10xUMIs.py!'''
def map_reads (inFile,ref_fasta,path):
    #sys.stderr.write('Mapping reads to %s' %(ref_fasta))
    out= path + '/allfinalreads.sam'
    os.system('%s --secondary=no -ax splice \
              %s %s > %s 2> ./minimap2_messages.txt'
              % (minimap2,ref_fasta, inFile, out))
    return out

'''Convert sam alignment file to bam'''
def sam_con (inFile):
    #sys.stderr.write('Converting sam to bam file')
    outname=inFile.split('.')[0]
    os.system('samtools view -S -b %s > %s.bam' %(inFile, outname))
    os.system('samtools sort -o %s.sorted.bam %s.bam' %(outname, outname))
    os.system('samtools index -b %s.sorted.bam' %(outname))


'''Assign mapped reads to gene names with featureCount--uses same parameters as Map10xUMIs.py!'''
def fcount (inFile, gtf):
    #sys.stderr.write('Getting gene assignments for mapped reads')
    os.system('featureCounts --primary -R CORE -Q 20 -L -t exon -g gene_id -a %s -o finalcountmtx_featureCount.txt %s' %(gtf, inFile))
#-L denotes long-reads, -Q 20 only assigns alignments with a MAPQ score greater than 20


def main(input_path, path, ref_genome, gtf_file):
    catfile=concat_fasta(input_path, path)
    #sys.stderr.write('Trimming and concatenating fastas')
    samfile=map_reads(catfile, ref_genome, path)
    sam_con(samfile)
    bamfile=path+'/'+'allfinalreads.bam' 
    fcount(bamfile, gtf)

main(input_path, path, ref_genome, gtf)
