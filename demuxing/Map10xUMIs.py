#!/usr/bin/env python3

import argparse
import editdistance
import numpy as np
import os
import sys
import mappy as mm
import shutil 
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path to individual cell .fasta and _subs.fastq output files from Vollmers demultiplex_kmer script
                                                    #final.fasta files with mapped and collapsed UMIs will go into this folder as well
parser.add_argument('-o', '--output_path', type=str) #where you want temp output files and assignedreads directory to go 
parser.add_argument('-g', '--gtf_file', type=str) #path to reference genome gtf file for FeatureCount gene id assignments
parser.add_argument('-c', '--config_file', type=str) #path to config file for minimap2, racon, FeatureCount, and samtools if not in $PATH
parser.add_argument('-r', '--ref_genome', type=str) #reference GENOME for splice-aware minimap2 alignment


args = parser.parse_args()
#input_file=args.input_file
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
        
    possible = set(['minimap2', 'conk', 'consensus', 'racon', 'samtools','featureCounts']) #don't think I actually need psl2pslx and emtrey...
    inConfig = set() 
    for key in progs.keys(): 
        inConfig.add(key)
        if key not in possible: 
            raise Exception('Check config file')

    for missing in possible-inConfig: 
        if missing == 'consensus': 
            path = 'consensus.py' #change this to where 'consensus.py is located in your computer' 
        else:
            path = missing  
        progs[missing] = path 
        #sys.stderr.write('Using ' + str(missing)
                         #+ ' from your path, not the config file.\n')
    return progs

progs = configReader(config_file)
#poa = progs['poa'] 
minimap2 = progs['minimap2']
racon = progs['racon']
#consensus = progs['consensus']
samtools = progs['samtools']
featureCounts = progs['featureCounts']


'''Parse each cell fasta, trim barcode and UMI sequences from each read and concat into one fasta file'''
def trim (inFile):
    readdict={}
    filename=inFile.split('/')[-1]
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
        if file.startswith('cell_') and file.endswith('.fasta'):
            fileList.append(file)       #fileList is a list containing files with '.fasta' extension in the 'input_path' directory in arbitrary order
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])):   #sorted files in increasing number based on cell_number
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
def fcount (inFile, gtf, path):
    #sys.stderr.write('Getting gene assignments for mapped reads')
    outdir = path + '/'
    os.system('featureCounts --primary -R CORE -Q 20 -L -t exon -g gene_id -a %s -o %smysample_featureCount.txt %s' %(gtf, outdir, inFile))
#-L denotes long-reads, -Q 20 only assigns alignments with a MAPQ score greater than 20

'''Group assigned reads by cell barcode and gene assignment'''
def group_reads (inFile, input_path, output_path):
    if not os. path. isdir(path + '/' + 'assignedreads'):
        #shutil.rmtree(path + '/' + 'assignedreads')
    #else:
        os.mkdir(path + '/' + 'assignedreads')
 
    previous_cell=''
    dict={}
    ass_dict={} #haha
    fileList=[]
    bcdict={}
    fullnamedict={}
    for line in open(inFile): #inFile is allreads.bam.featurecounts
        name=line.rstrip().split('\t')[0] #full name containing readname(header) and cell_#_16bpBarcode
        rname=name.split('_')[0]          #readname(header)
        if line.rstrip().split('\t')[3]!='NA'and rname not in ass_dict: #just take the primary assignment
            ass_dict[rname]=line.rstrip().split('\t')[3] #assignment dictionary at that fullname stores the gene id for that read--reads are being duplicated! Need to figure out how to fix this!
        cell=name.split('_')[2] #takes cell number from each read name (doesn't care about assignment)
        bc=name.split('_')[3] #takes the barcode from each read name (doesn't care about assignment)
        bcdict[cell]=bc #stores barcode information/sequence (dict_value) under cell number(dict_key) in the barcode dictionary
        #final = open(path + '/' + cell + '.new.fasta', 'w')
        if int(line.rstrip().split('\t')[2])==1:
            if cell!=previous_cell:   #check if the same cell 
                previous_cell=cell    #previous_cell is the cell number
                dict[cell]=[]
                bcdict[cell]=bc
                if name not in dict.get(cell): #my clumsy way of preventing read duplicates and taking the first gene_id in featurecounts (could be better)
                    dict[cell].append(name)
            else:
                if name not in dict.get(cell):
                    dict[cell].append(name)
                    #dict stores the cell number (dict_key) and then all the reads that had gene assignments (in a list) (dict_value)
    dictlength=int(cell) #dictionary has as many entries as the last cell counted that had gene assignments (some cells don't have gene assignments)
    for i in range(dictlength+1): #line160-163 is for cells that do not have gene assignments
        i=str(i)
        if i not in dict:
            dict[i]=[] #for cells without gene assignments make those cells empty lists in the dict
    for i in range(len(dict)):
        readlist=[]
        i=str(i)
        barcode=bcdict[i] #attach barcode info to cell 
        if os.path.exists(path + '/assignedreads' + '/'+'cell_'+ i + '_'+barcode+'.new.fasta'):
            os.remove(path + '/assignedreads' + '/'+'cell_'+ i + '_'+barcode+'.new.fasta')
        final = open(path + '/assignedreads' + '/'+'cell_'+ i + '_'+barcode+'.new.fasta', 'w') #open file with cell number and barcode
        file=input_path+'/'+'cell_'+i+'_'+barcode+'.fasta'
        seqdict={} 
        for line in open(file): #takes the reads from each cell fasta
            if line.startswith('>'):
                #line=line.rstrip()
                rootname=line.split('_')[0][1:]
                fullname=line.rstrip()[1:] #fetches the full name for each read in these files (so we can determine best consensus read for polishing)
                seqdict[rootname]=[]
                fullnamedict[rootname]=fullname #fullnamedict is cumulative across all single cell fastas 
            else:
                sequence=line.rstrip()
                seqdict[rootname]=sequence #fetches sequence for each read in this file
        for x in dict[i]:
            x=x.split('_')[0]
            name=fullnamedict[x]
            seq=seqdict[x]
            ass=ass_dict[x]
            final.write('%s\t%s\t%s\n' % (name,ass,seq))
        final.close()
    for file in os.listdir(path+'/assignedreads'):
        if file.endswith('new.fasta'):
            fileList.append(file)
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])): #creates cell_#_BC.fasta.new.sorted file with reads sorted by gene_id
        with open(path+'/assignedreads/'+file)as open_file:
            if os.path.exists(path + '/assignedreads/'+ file +'.sorted'):
                os.remove(path + '/assignedreads/'+ file +'.sorted')
            sorted_file=open(path + '/assignedreads/'+ file +'.sorted', 'w')
            rows=open_file.readlines()
            sorted_rows = sorted(rows, key=lambda x: int(x.split()[1][7:]), reverse=False)
            for row in sorted_rows:
                sorted_file.write(row)

'''Evaluate UMIs for reads that have same barcode and cell ID, if UMIs are L=<1, polish using the 'best' read as the template'''
def eval_UMIs(path, input_path):
    fileList=[] 
    for file in os.listdir(path+'/assignedreads/'):
        if file.endswith('.new.fasta.sorted'):
            fileList.append(file)
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])): 
        cell=file.split('_')[1]  #cell is the cell number
        barcode=(file.split('_')[2]).split('.')[0]   #barcode is the 16bp barcode sequence
        if os.path.exists(input_path + '/'+'cell_'+ cell + '_'+barcode+'.final.fasta'):
            os.remove(input_path + '/'+'cell_'+ cell + '_'+barcode+'.final.fasta')
        finalfasta = open(input_path + '/'+'cell_'+ cell + '_'+barcode+'.final.fasta', 'a') #if you already have this file it will be APPENDED so if you are doing this more than once be mindful
        genedict={}
        previous_gene=''

        for line in open(path+'/assignedreads/'+file):
            line=line.rstrip()
            gene=line.split('\t')[1]
            name=line.split('\t')[0]
            sequence=line.split('\t')[2]
            UMI=sequence[len(sequence)-17::-1][0:12]
            if gene != previous_gene:
                previous_gene=gene
                genedict[gene]=[]
                genedict[gene].append(name)
                genedict[gene].append(sequence)
                genedict[gene].append(UMI)
            else:
                genedict[gene].append(name)
                genedict[gene].append(sequence)
                genedict[gene].append(UMI) #creates genedict that has gene_id keys storing the names, seqs, and UMIs of associated reads (in that order)
        #genedict{} stores contents of 'new.fasta.sorted'; key is 'gene_id', values as ['readname1','seq1','UMI1'(,'readname2','seq2','UMI2'...)]
        for key,value in genedict.items(): #key is referring to gene id
            length=len(list(filter(None, value)))
            if length < 6: #if there is only one read for this gene in the cell...
                finalfasta.write('>%s\n%s\n' % (genedict[key][0], genedict[key][1]))

            elif length >= 6: #if gene has two or more reads associated with it, check for matching UMIs between them
                group={}
                y=0
                id=-1
                for i in range(2,length,3): #creates little readname subgroups with matched UMIs 
                    for j in range(i+3,length,3):
                        id+=1
                        group[id]=[]
                        if editdistance.eval(genedict[key][i],genedict[key][j]) <=1:
                            condition='no' 
                            while genedict[key][i-2] and genedict[key][j-2] not in group[y]: #as groupdict grows, check each previous entry for a match
                                condition='no'
                                if genedict[key][i-2] in group[y]:
                                    group[y].append(genedict[key][i-2])
                                    group[y].append(genedict[key][i-1])
                                    group[y].append(genedict[key][j-2])
                                    group[y].append(genedict[key][j-1])
                                    condition='yes'
                                    break #if match is found, append the values and break loop
                            
                                elif y>=id: #once the length of the dictionary has been reached, break the loop
                                    break

                                y += 1 #loop through the dictionary until either of the above conditions are met

                                if condition == 'no':
                                    group[id]=[genedict[key][i-2], genedict[key][i-1],genedict[key][j-2],genedict[key][j-1]]
                                    #if no matches found after scanning all previous entries, add new entry with the values
 
                for entry in group: #clean up group dictionary and consolidate entries
                    for z in range(entry+1,len(group),1):
                        for a in range(0,len(group[entry]),1):
                            if group[entry][a] in group[z]:
                                group[z]=[]
                #print(group)

                for entry in group:
                    if group[entry]:
                        refgroup={}
                        leng=len(group[entry])
                        for i in range(0,leng,2):
                            name=group[entry][i]
                            seq=group[entry][i+1]
                            refgroup[name]=seq #creates readname:seq dictionary for all matched reads
                        readname,readseq=make_consensus(refgroup,path, input_path, cell, barcode) #call make_consensus function--defined below
                        finalfasta.write('%s\n%s\n' % (readname, readseq)) #writes consensus read to out file
                matchlist=[]
                for entry in group:
                    if group[entry]:
                        for c in range(0,len(group[entry]),1):
                            matchlist.append(group[entry][c])

                for b in range(0,length,3): #writes reads under the gene_id that had NO UMI matches to final cell fasta
                    if genedict[key][b] not in matchlist:
                        finalfasta.write('>%s\n%s\n' % (genedict[key][b], genedict[key][b+1]))
        finalfasta.close()

'''Determines which read in the group has most subreads, aligns all fastqs to that read with minimap and polishes with racon'''
def make_consensus(group, path, input_path, cell, barcode):
    previoussub=0
    if os.path.exists('targettemp.fasta'): #make temp target fasta for each batch of reads to be aligned and polished
        os.remove('targettemp.fasta')
    if os.path.exists('querytemp.fastq'):
        os.remove('querytemp.fastq')
    targetFa = open(path+'/targettemp.fasta', 'w') #make temp fasta for consensus fasta with most subreads, will be contig for align and correct 
    queryFq = open(path+'/querytemp.fastq','a') #make temp fastq file for query sequences
    for key in group:
        subreads=int(key.split('_')[3])
        if subreads > previoussub:
            previoussub=subreads
            best=key
    targetFa.write('>%s\n%s\n' % (best,group[best]))
    targetFa.close()
    add=-1
    linenum=0
    for line in open(input_path+'/'+'cell_'+str(cell)+'_'+barcode+'_subs.fastq'): #go fetch the fastqs for each consensus read in the group
        linenum +=1
        if line.startswith('@') and linenum%4==1:
            name=line.rstrip()[1:]
        if linenum%4==2:
            seq=line.rstrip()
        if linenum%4==3:
            strand=line.rstrip()
        if linenum%4==0:
            qual=line.rstrip()

            for key in group:
                if name in key:
                    add +=1
                    #print(name,seq,strand,qual)
                    queryFq.write('@%s_%s\n%s\n%s\n%s\n' %(name,str(add),seq,strand,qual))
    queryFq.close()

    #print(barcode)
    #print(cell)
    refFasta=path+'/targettemp.fasta'
    inFastq=path+'/querytemp.fastq'
    overlap=path+'/tempalignment.sam'
    finalout=path+'/temppolish.fasta'

    os.system('%s --secondary=no -ax map-ont %s %s > %s 2> ./minimap2_messages.txt' % (minimap2, refFasta, inFastq, overlap)) #align to best read

    os.system('%s -q 5 -t 1 %s %s %s > %s 2>./racon_messages.txt' % (racon, inFastq, overlap, refFasta, finalout))
    #uses best consensus sequence as contig against which to align and error-correct all of the fastqs

    for line in open(path+'/temppolish.fasta'):
        line=line.rstrip()
        if line.startswith('>'):
            readname=line.split(' ')[0]
        else:
            readseq=line

    return readname,readseq


def main(input_path, path, ref_genome, gtf_file):
    catfile=concat_fasta(input_path, path)
    #sys.stderr.write('Trimming and concatenating fastas')
    samfile=map_reads(catfile, ref_genome, path)
    sam_con(samfile)
    bamfile=path+'/'+'allreads.bam' 
    fcount(bamfile, gtf,path)
    countfile=path+'/'+'allreads.bam.featureCounts'
    group_reads(countfile,input_path, path)
    eval_UMIs(path,input_path)


main(input_path, path, ref_genome, gtf)
