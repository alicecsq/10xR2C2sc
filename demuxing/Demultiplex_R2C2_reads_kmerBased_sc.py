#!/usr/bin/env python3
'''Sheridan 2021 edits'''

import sys
import argparse
import editdistance as ld #this is for L distance
from tqdm import tqdm
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_fasta_file', type=str)
parser.add_argument('-n', '--tenx_index_file', type=str)
parser.add_argument('-o', '--output_path', type=str)

args = parser.parse_args()
input_file = args.input_fasta_file
bc_file = args.tenx_index_file
output_path = args.output_path + '/'

def read_fasta(inFile, bc): #reads input fasta and builds readDict if not barcode file (bc=False) or bcKmerDict if barcode file indicated (bc=True)
    readDict = {}
    bcKmerDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            readDict[line[1:]] = '' 
            lastHead = line[1:] 
        else: 
            readDict[lastHead] += line 
            if bc: #put a condition in the read_fasta file that denotes whether its the barcode file or not
                kmers = get_kmer_list(line) #kmers is the list containing 4mers and positions
                for kmer in kmers: 
                    if not bcKmerDict.get(kmer): #if kmer is not in bckmerDict
                        bcKmerDict[kmer] = [(lastHead, readDict[lastHead])] #bcKmerDict at this kmer stores the list of the barcode name and sequence
                    else:
                        bcKmerDict[kmer].append((lastHead, readDict[lastHead])) #theoretically, kmers CAN be shared between barcodes because some of these 10x barcodes are only L=1 
    if bc:
        return bcKmerDict
    else:
        return readDict

def get_kmer_list(bc): #operated on each sequence in the freqbarcodes file
    kmer_list = [bc[:4] + '___', #generates a list of four variables that encode 4mer sequence and position with '_' characters
                 '_' + bc[4:8] + '__',
                 '__' + bc[8:12] + '_',
                 '___' + bc[12:16]]
    return kmer_list

def demultiplex(reads, bcKmerDict): #reads=readname:barcode dictionary, bcKmerDict=kmerkey:(barcode name,barcode sequence) dictionary
    indexed_reads, counter = {}, 0
    totalcount=0
    for read, complete_sequence in tqdm(reads.items()): #reads.items() separates the reads dictionary into key:value pairs, so the len of reads.items() is the number of total reads
        # if counter % 10000 == 0:
        #     print(str(counter) + ' of ' + str(len(reads)))
        sequence = complete_sequence[:16] 
        totalcount += 1
        bc_set = set()
        readKmers = get_kmer_list(sequence) #for each read, generate kmer keys from the barcode associated with the read (again this is the same matching strategy that was used for the splint UMI kmers)
        for kmer in readKmers:
            if bcKmerDict.get(kmer): #if one of the kmer keys from the read matches one from the 1500 most freq dictionary...
                barcodes = bcKmerDict[kmer]
                for b in barcodes: 
                    bc_set.add(b)
        dists, bc_match, maxDist = [], '', 3
        for barcode in bc_set:
            dist = ld.eval(barcode[1], sequence)
            dists.append((barcode, dist))
        if dists:
            sorted_bc_list = sorted(dists, key=lambda x: x[1]) #sorting matched barcodes from smallest L distance to largest

            if sorted_bc_list[0][1] < maxDist and sorted_bc_list[0][1] < sorted_bc_list[1][1]-1: #and is less than the L distance of the next smallest match (i.e. is it a UNIQUE match) then...
                bc_match = sorted_bc_list[0][0][0]
                counter += 1 
            read += '|' + bc_match #read becomes read name plus a pipe and the name of the matched barcode
            indexed_reads[read] = complete_sequence #indexed reads at key 'readname|matched_bc' holds the value of the read sequence
    percentdemux=(counter/totalcount)*100
    print('Pipeline demuxed %s percent of %s total reads' %(percentdemux,totalcount)) 
    return indexed_reads

def write_fasta_file(path, reads):
    out = open(path + '/allbckmer_demuxed.fasta', 'w+')
    for name, sequence in reads.items():
        out.write('>%s\n%s\n' %(name, sequence))


'''def makeHist(counts):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mplpatches
    # plt.style.use('BME163')

    plt.figure(figsize = (6, 2))
    hist = plt.axes([0.1, 0.1, 8/10, 4/5], frameon = True)

    for i in range(0,len(counts)):
        width, height = 1, counts[i]
        bar = mplpatches.Rectangle((i, 0), width, height,
                                   lw=0, fc='black')
        hist.add_patch(bar)

    hist.set_xlim(0, 1800)
    hist.set_ylim(0, 1500)
    hist.set_title('Nanopore barcode counts')
    hist.set_xlabel('Barcode')
    hist.set_ylabel('Count')
    hist.tick_params(axis='both',which='both',\
                     bottom='on', labelbottom='on',\
                     left='on', labelleft='on',\
                     right='off', labelright='off',\
                     top='off', labeltop='off')
    plt.savefig('combined_ro_hist.png', bbox_inches="tight")'''



bcKmerDict = read_fasta(bc_file, True)
reads = read_fasta(input_file, False)
indexed_reads=demultiplex(reads, bcKmerDict)
write_fasta_file(output_path, indexed_reads)

#histogram of demuxed and non demuxed reads, subreads per raw read
matchedsubreadlist=[]
subreadlist=[]
for item in indexed_reads.items():
    name=item[0].split('_')
    barcodevet=item[0].split('|')
    if barcodevet[2]=='':
        subreadlist.append(int(name[3]))
    else:
        matchedsubreadlist.append(int(name[3]))
matchedsubreadlist=sorted(matchedsubreadlist, reverse=True)
subreadlist=sorted(subreadlist, reverse=True)


import matplotlib.pyplot as plt
plt.hist(subreadlist, bins=100, alpha=0.5, label="non-demuxed")
plt.hist(matchedsubreadlist, bins=100, alpha=0.5, label="demuxed")
plt.legend(loc='upper right')
plt.yscale('log', nonposy='clip')
plt.xlim(-10,200)
plt.xlabel("Number of Subreads", size=14)
plt.ylabel("Read Count", size=14)
plt.title("Subreads per Raw Read Distribution, BC=20k")
plt.savefig('subreaddistallbc.png', bbox_inches="tight")

unmuxmedian=np.median(subreadlist)
demuxmedian=np.median(matchedsubreadlist)
print('The median subread number for demuxed reads is %s, the median subread number for non-demuxed reads is %s.' %(demuxmedian,unmuxmedian))


'''n = 0
toptotal=0
for count in counts[:1500]:
    print('>barcode_' + str(n) + '_' + str(count[0]))
    print(count[1])
    toptotal += count[0]
    n += 1
makeHist(counts)'''



#to look at percent of merged_UMI reads that are demultiplexed
totalread=0
matches=0
for item in indexed_reads.items():
    barcodevet=item[0].split('|')
    name=item[0].split('|')
    beef=name[0].split('_')
    if len(beef)>5:
        totalread += 1
    if len(beef)>5 and barcodevet[2] != '':
        matches += 1
UMImatched=(matches/totalread)*100
print('out of %s splint-UMI-merged reads reads, %s percent are demultiplexed' %(totalread,UMImatched))



#un-annotate to get a count of the plus vs minus strands that were demuxed
minusreads={}
minustotalcount = 0
minusemptycount = 0
for item in indexed_reads.items():
    if 'minus' in item[1]:
        minusreads[item[0]]=item[1]
        minustotalcount += 1
for name in minusreads:
    jeff=name.split('|')
    if jeff[2]=='':
        minusemptycount += 1
minusdemux=100-((minusemptycount/minustotalcount)*100)
print('Out of %s minus-strand reads, %s percent were demultiplexed.' %(minustotalcount,minusdemux))

plusreads={}
plustotalcount = 0
plusemptycount = 0
for item in indexed_reads.items():
    if 'plus' in item[1]:
        plusreads[item[0]]=item[1]
        plustotalcount += 1
for name in plusreads:
    suzie=name.split('|')
    if suzie[2]=='':
        plusemptycount += 1
plusdemux=100-((plusemptycount/plustotalcount)*100)
print('Out of %s plus-strand reads, %s percent were demultiplexed.' %(plustotalcount,plusdemux))
