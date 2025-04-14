#!/usr/bin/env python

import allel
import numpy as np

import re
import sys
import argparse 
import pandas as pd

#from plotnine import *


# vcf='pileup_all_vars.vcf.gz'
# ref="reference.fasta"
# #gff="reference.gff"
# rangesbed="TBscheme.insert.bed"
# depthfile="all_depth.txt.gz"
# outfile="TB001_amplicon"


parser = argparse.ArgumentParser()
parser.add_argument('--fa', '-f', help='reference seq (fasta)')
parser.add_argument('--bed', '-b', help='primer / amplicon positions (bed)')
parser.add_argument('--vcf', '-v', help='variant calls (vcf)')
parser.add_argument('--out', '-o', help='outfile prefix')

args = parser.parse_args()
fafile = args.fa
bedfile = args.bed
vcffile = args.vcf
outfile = args.out


ranges=pd.read_csv(bedfile, sep='\t',header=None)
callset = allel.read_vcf(vcffile, fields='*')

gt = allel.HaplotypeArray(callset['calldata/GT'][:,:,0])
ac = gt.count_alleles()
pos = callset['variants/POS']
chrom = callset['variants/CHROM'].tolist()[0]
samples = callset['samples']

nranges=len(ranges)
nsamples=len(samples)

pi=pd.DataFrame([[-1.0]*(nsamples+2)]*nranges,columns=['start', 'end']+samples.tolist())

for i in range(0,len(ranges)):
    st = ranges[1][i]
    en = ranges[2][i]
    len=en-st
    inpos = np.logical_and((pos >= st),(pos <= en))

    if sum(inpos)>0:
        ingts = gt[inpos,:]
        locpi = (ingts[0] / len).round(4).tolist()
        #print([st,en]+locdepth+locpi+loccov)
    else:
        locpi = [-1.0]*nsamples
    pi.iloc[i,:] = [st,en]+locpi

pi.to_csv(outfile+"-pi.tsv",sep="\t")

pim = pd.melt(pi,id_vars=["start","end"],value_name="pi",var_name="sample")




#piplot = ggplot(pim, aes(x="sample",y="pi"))+ geom_violin() +coord_flip()
#piplot.save((outfile+"_pi.png"),width=150,height=150,units="mm")
