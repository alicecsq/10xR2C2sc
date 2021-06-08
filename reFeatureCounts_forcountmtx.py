#!/usr/bin/env python3
'''Script for re-mapping individual cell reads for FeatureCounts and Seurat formatting'''
'''dependencies:minimap2'''



import argparse
import editdistance
import numpy as np
import os
import sys
import mappy as mm
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path to individual cell .final.fasta UMI-merged/polished files
parser.add_argument('-o', '--output_path', type=str) #where you want 'cellsams' directory to go 
parser.add_argument('-c', '--config_file', type=str) #config file containing paths to executables if not in $PATH
parser.add_argument('-r', '--ref_genome', type=str) #reference GENOME for minimap2 alignment (same used for Map10xUMIs.py)


args = parser.parse_args()
input_path=args.input_path
path=args.output_path
config_file=args.config_file
ref_genome=args.ref_genome

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


'''Read each cell fasta in input path, create temp trimmed fasta for minimap2 alignment'''

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


def map_reads(name,inFile,ref_genome,path):
    #sys.stderr.write('Mapping reads to %s' %(ref_fasta))
    out= path + '/cellsams/'+ name +'.sam'
    os.system('%s --secondary=no -ax splice %s %s > %s' % (minimap2, ref_genome, inFile, out))
    #return out

def concat_fasta(input_path, path, ref_genome):
    #sys.stderr.write('Concatenating reads')
    fileList=[]
    for file in os.listdir(input_path):
        #final=open(path+'/temptrimmed.fasta','w')
        if '.final.fasta' in file:
            fileList.append(file)
    #print(sorted(fileList, key=lambda x: int(x.split('_')[1])))
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])):
        name='cell'+ '_'+str((file.split('.')[0]).split('_')[1])
        final=open(path+'/temptrimmed.fasta','w')
        temp=path+'/temptrimmed.fasta'
        fasta_file = input_path + '/' + file
        readdict=trim(fasta_file)
        for i in readdict:
            sequence=readdict[i] 
            final.write('%s\n%s\n' % (i,sequence))
        #final.close()
        map_reads(name,temp,ref_genome,path)
    #return path+'/allfinalreads.fasta'


def main(input_path, path, ref_genome):
    if not os. path. isdir(path+'cellsams'):
        #shutil.rmtree(path+'assignedreads')
    #else:
        os.mkdir(path+'cellsams')

    catfile=concat_fasta(input_path, path, ref_genome)

main(input_path, path, ref_genome)
