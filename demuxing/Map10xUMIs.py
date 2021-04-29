#!/usr/bin/env python3

import argparse
import editdistance
import numpy as np
import os
import sys
import mappy as mm 
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str)
parser.add_argument('-i', '--input_file', type=str)
parser.add_argument('-o', '--output_path', type=str)
#parser.add_argument('-g', '--gtf_file', type=str)
#parser.add_argument('-c', '--config_file', type=str)
#parser.add_argument('-r', '--ref_genome', type=str)


args = parser.parse_args()
input_file=args.input_file
input_path=args.input_path
path=args.output_path
#config_file=args.config_file
#ref_genome=args.ref_genome
#gtf=args.gtf_file


'''Parse config'''
def configReader(configIn):
    progs = {}
    for line in open(configIn): 
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t') 
        progs[line[0]] = line[1] 
        
    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'samtools','featureCount', 'psl2pslx']) #don't think I actually need psl2pslx and emtrey...
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

'''progs = configReader(config_file)
poa = progs['poa'] 
minimap2 = progs['minimap2']
racon = progs['racon']
consensus = progs['consensus']
samtools = progs['samtools']
featureCount = progs['featureCount']'''


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
    final = open(path + '/allreads.fasta', 'w')
    fileList=[]
    for file in os.listdir(input_path):
        if '.fasta' in file:
            fileList.append(file)
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])):
        fasta_file = input_path + '/' + file
        readdict=trim(fasta_file)
        for i in readdict:
            sequence=readdict[i] 
            final.write('%s\n%s\n' % (i,sequence))
    return path+'/allreads.fasta'


'''Map reads in the concatenated fasta file with minimap2'''
def map_reads (inFile,ref_fasta,path):
    #sys.stderr.write('Mapping reads to %s' %(ref_fasta))
    out= path + '/allreads.sam'
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


'''Assign mapped reads to gene names with featureCount'''
def fcount (inFile, gtf):
    #sys.stderr.write('Getting gene assignments for mapped reads')
    os.system('featureCounts -R CORE -Q 20 -L -t exon -g gene_id -a %s -o mysample_featureCount.txt %s' %(gtf, inFile))
#-L denotes long-reads, -Q 20 only assigns alignments with a MAPQ score greater than 20

'''Group assigned reads by cell barcode and gene assignment'''
def group_reads (inFile, input_path, output_path):
    if not os. path. isdir(path+'assignedreads'):
        os.mkdir(path+'assignedreads')
    previous_cell=''
    dict={}
    ass_dict={} #haha
    fileList=[]
    bcdict={}
    fullnamedict={}
    for line in open(inFile):
        name=line.rstrip().split('\t')[0]
        if line.rstrip().split('\t')[3]!='NA':
            ass_dict[name]=line.rstrip().split('\t')[3]
        cell=name.split('_')[2]
        bc=name.split('_')[3]
        bcdict[cell]=bc
        #final = open(path + '/' + cell + '.new.fasta', 'w')
        if int(line.rstrip().split('\t')[2])==1:
            if cell!=previous_cell:
                previous_cell=cell
                dict[cell]=[]
                bcdict[cell]=bc
            else:
                dict[cell].append(name)
    dictlength=int(cell)
    for i in range(dictlength+1):
        i=str(i)
        if i not in dict:
            dict[i]=[]
    for i in range(len(dict)):
        i=str(i)
        barcode=bcdict[i]
        final = open(path + '/assignedreads' + '/'+'cell_'+ i + '_'+barcode+'.new.fasta', 'w')
        file=input_path+'/'+'cell_'+i+'_'+barcode+'.fasta'
        seqdict={}
        for line in open(file):
            if line.startswith('>'):
                rootname=line.split('_')[0][1:]
                fullname=line.rstrip()[1:]
                seqdict[rootname]=[]
                fullnamedict[rootname]=fullname
            else:
                sequence=line.rstrip()
                seqdict[rootname]=sequence
        for x in dict[i]:
            x=x.split('_')[0][1:]
            print(x)
            #if x in fullnamedict:
                #print(x, fullnamedict[x])
                #name=fullnamedict[x]
                #seq=seqdict[name]
                #ass=ass_dict[x]
            #final.write('%s\n%s\n%s\n' % (name,seq,ass))
    #print(sorted(fileList, key=lambda x: int(x.split('_')[1])))

        #for x in dict[i]:
            #print(x,ass_dict[x]) 
 
    '''for file in os.listdir(input_path):
        if '.fasta' in file:
            fileList.append(file)
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])):
        fasta_file = input_path + '/' + file
        readdict=trim(fasta_file)
        for i in readdict:
            sequence=readdict[i]
            final.write('%s\n%s\n' % (i,sequence))'''

 
    '''final_UMI_only = open(path + '/' + process_count + '.UMI_only.fasta', 'w')
    final_subreads = open(path + '/' + process_count + '.merged.subreads.fastq', 'w')
    matched_reads = open(path + '/' + process_count + '.matched_reads.txt', 'w')''' 


'''def main(input_path, path, ref_genome, gtf_file):
    catfile=concat_fasta(input_path, path)
    #sys.stderr.write('Trimming and concatenating fastas')
    samfile=map_reads(catfile, ref_genome, path)
    sam_con(samfile)
    bamfile=path+'/'+'allreads.bam' 
    fcount(bamfile, gtf)
    countfile=path+'/'+'allreads.bam.featureCounts'
    group_reads(countfile,input_path, path)

main(input_path, path, ref_genome, gtf)'''

group_reads(input_file,input_path, path)



#catfile=concat_fasta(input_path,path)
#map_reads(catfile, ref_genome, path)

#for line in open(catfile):
    #print(line.strip())
