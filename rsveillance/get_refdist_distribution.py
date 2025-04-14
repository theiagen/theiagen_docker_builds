#!/usr/bin/env python

import re
import sys
import argparse 
from math import floor, ceil 

from statistics import mean
import allel
import numpy as np
import pandas as pd

# vcf='pileup_all_vars.vcf.gz'
# ref="reference.fasta"
# #gff="reference.gff"
# rangesbed="TBscheme.insert.bed"
# depthfile="all_depth.txt.gz"
# outfile="TB001_amplicon"

parser = argparse.ArgumentParser()
parser.add_argument('--window', '-w', help='depth window size',default=1000)
parser.add_argument('--bed', '-b', help='primer / amplicon positions (bed)')
parser.add_argument('--gff', '-g', help='gene positions (gff)')
parser.add_argument('--vcf', '-v', help='variant calls (vcf)')
parser.add_argument('--out', '-o', help='outfile prefix')
parser.add_argument('--sample', '-s', help='outfile prefix')

args = parser.parse_args()
winsize = int(args.window)
vcffile = args.vcf
bedfile = args.bed
gfffile = args.gff
sample = args.sample
outfilename = args.out


outfile = open(outfilename,'w')

#get ranges from bed / gff file
ranges=list()
if bedfile is not None:
    rangetab=pd.read_csv(bedfile, sep='\t')
    for i in range(0,len(rangetab)):
        st = int(rangetab.iloc[i,1])
        en = int(rangetab.iloc[i,2])
        name = rangetab.iloc[i,3]
        ranges.append((st,en,name))

elif gfffile is not None:
    getname = re.compile(r'.*?gene=([^;]*);.*')

    rangetab=pd.read_csv(gfffile, sep='\t', comment="#")
    for i in range(0,len(rangetab)):
        st = int(rangetab.iloc[i,3])
        en = int(rangetab.iloc[i,4])
        type = rangetab.iloc[i,2]
        namestr = rangetab.iloc[i,8]
        if type=="CDS":
            search = getname.search(namestr)
            if search is not None:
                name = search.group(1)
            else:
                name = "NA"
            ranges.append((st,en,name))


#read in VCF 
callset = allel.read_vcf(vcffile, fields='*')
gt = allel.HaplotypeArray(callset['calldata/GT'][:,:,0])
ac = gt.count_alleles()
chrom = callset['variants/CHROM'].tolist()[0]
samples = callset['samples']
pos = callset['variants/POS']
nsamples=len(samples)
    
def getpiregion(st,en):
    inpos = np.logical_and((pos >= st),(pos <= en))
    len=(en-st)+1
    if sum(inpos)>0:
        ingts = gt[inpos,:]
        locpi = (np.sum(ingts,0) / len).round(4).tolist()
        #print([st,en]+locdepth+locpi+loccov)
    else:
        locpi = [-1.0]*nsamples
    return locpi    


#get wins from bed if bedfile given
if len(ranges) >0:
    nranges=len(ranges)
    for st,en,nm in ranges:
        locpi = getpiregion(st,en)
        #pistring = "\t".join(map(str,locpi))
        #print("{}\t{}\t{}".format(st,en,pistring), file=outfile)
        for spi in zip(samples,locpi):
            s,pi = spi
            print("{}\t{}\t{}\t{}\t{}".format(st,en,nm,s,pi), file=outfile)
            
else:
    nranges = ceil(max(pos)/winsize)
    for st in range(1,nranges*winsize,winsize):
        en = st+winsize-1
        locpi = getpiregion(st,en)
        #pistring = "\t".join(map(str,locpi))
        #print("{}\t{}\t{}".format(st,en,pistring), file=outfile)
        for spi in zip(samples,locpi):
            s,pi = spi
            print("{}\t{}\t{}\t{}".format(st,en,s,pi), file=outfile)



