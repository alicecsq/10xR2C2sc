#!/usr/bin/env python3

import os
import sys
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path to individual cell .fasta and _subs.fastq output files from Vollmers demultiplex_kmer script
                                                    #final.fasta files with mapped and collapsed UMIs will go into this folder as well
parser.add_argument('-o', '--output_path', type=str) #where you want temp output files and assignedreads directory to go
#parser.add_argument('-g', '--gtf_file', type=str) #path to reference genome gtf file for FeatureCount gene id assignments
#parser.add_argument('-c', '--config_file', type=str) #path to config file for minimap2, racon, FeatureCount, and samtools if not in $PATH
#parser.add_argument('-r', '--ref_genome', type=str) #reference GENOME for splice-aware minimap2 alignment


args = parser.parse_args()
#input_file=args.input_file
input_path=args.input_path
path=args.output_path
#config_file=args.config_file
#ref_genome=args.ref_genome
#gtf=args.gtf_file

def sort_files(input_path,path):
    fileList=[]
    for file in os.listdir(input_path):
        fileList.append(file)
        #cellnum=(file.split('_')[1]).split('.')[0]
        #print(name)
    zeroes=len(str(int(len(fileList)/10)))
    #print(zeroes)
    i=0
    for sam in sorted(fileList, key=lambda x: int((x.split('_')[1]).split('.')[0])):
        cellnum=(sam.split('_')[1]).split('.')[0]
        #print(cellnum)
        if i==len(fileList):
            break
        elif i < 10:
            shutil.copyfile(input_path+sam,path+'cell_'+('0'*zeroes)+cellnum+'.sam')
        elif i >=10 and i < 100:
            zeroes=(len(str(int(len(fileList)/10))))-1
            shutil.copyfile(input_path+sam,path+'cell_'+('0'*zeroes)+cellnum+'.sam')
        elif i >= 100 and i < 1000:
            zeroes=(len(str(int(len(fileList)/10))))-2
            shutil.copyfile(input_path+sam,path+'cell_'+('0'*zeroes)+cellnum+'.sam')
        elif i >= 1000 and i < 10000:
            zeroes=(len(str(int(len(fileList)/10))))-3
            shutil.copyfile(input_path+sam,path+'cell_'+('0'*zeroes)+cellnum+'.sam')
        elif i >= 10000 and i < 100000:
            zeroes=(len(str(int(len(fileList)/10))))-4
            shutil.copyfile(input_path+sam,path+'cell_'+('0'*zeroes)+cellnum+'.sam')
        i += 1 

if os. path. isdir(path+'/cellsamssorted'):
  os.mkdir(path+'/cellsamssorted')
sort_files(input_path,path)
