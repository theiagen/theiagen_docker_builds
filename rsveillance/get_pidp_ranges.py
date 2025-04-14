#!/usr/bin/env python

import allel
import numpy as np

import re
import sys
import argparse 
import pandas as pd

from plotnine import *


# vcf='pileup_all_vars.vcf.gz'
# ref="reference.fasta"
# #gff="reference.gff"
# rangesbed="TBscheme.insert.bed"
# depthfile="all_depth.txt.gz"
# outfile="TB001_amplicon"


parser = argparse.ArgumentParser()
parser.add_argument('--fa', '-f', help='reference seq (fasta)')
#parser.add_argument('--gff', '-g', help='primer scheme (bed)')
parser.add_argument('--bed', '-b', help='primer / amplicon positions (bed)')
parser.add_argument('--vcf', '-v', help='variant calls (vcf)')
parser.add_argument('--depth', '-d', help='depth file (samtools / txt)')
parser.add_argument('--out', '-o', help='outfile prefix')

args = parser.parse_args()
fafile = args.fa
bedfile = args.bed
vcffile = args.vcf
depthfile = args.depth
outfile = args.out


depth = pd.read_csv(depthfile, sep='\t')
ranges=pd.read_csv(bedfile, sep='\t',header=None)
callset = allel.read_vcf(vcffile, fields='*')

gt = allel.HaplotypeArray(callset['calldata/GT'][:,:,0])
ac = gt.count_alleles()
pos = callset['variants/POS']
chrom = callset['variants/CHROM'].tolist()[0]
samples = callset['samples']

nranges=len(ranges)
nsamples=len(samples)

covpc=pd.DataFrame([[-1]*(nsamples+2)]*nranges,columns=['start', 'end']+samples.tolist())
meandepth=pd.DataFrame([[-1]*(nsamples+2)]*nranges,columns=['start', 'end']+samples.tolist())
pi=pd.DataFrame([[-1.0]*(nsamples+2)]*nranges,columns=['start', 'end']+samples.tolist())

for i in range(0,len(ranges)):
    st = ranges[1][i]
    en = ranges[2][i]
    inpos = np.logical_and((pos >= st),(pos <= en))
    indepth = np.logical_and((depth['POS'] >= st),(depth['POS'] <= en))

    if sum(inpos)>0:
        ingts = gt[inpos,:]
        covbases = (depth[indepth].iloc[:,2:]>0).sum()
        locdepth = depth[indepth].iloc[:,2:].mean().tolist()
        locpi = (ingts[0] / covbases).tolist()
        loccov = (covbases / sum(indepth)).tolist()
        #print([st,en]+locdepth+locpi+loccov)
    else:
        locdepth = [0.0]*nsamples
        locpi = [-1.0]*nsamples
        loccov = [0.0]*nsamples
    pi.iloc[i,:] = [st,en]+locpi
    covpc.iloc[i,:] = [st,en]+loccov
    meandepth.iloc[i,:] = [st,en]+locdepth


meandepth.to_csv(outfile+"_meandp.tsv",sep="\t")
covpc.to_csv(outfile+"_covpc.tsv",sep="\t")
pi.to_csv(outfile+"_pi.tsv",sep="\t")


meandepthm = pd.melt(meandepth,id_vars=["start","end"],value_name="depth",var_name="sample")
covpcm = pd.melt(covpc,id_vars=["start","end"],value_name="covpc",var_name="sample")
pim = pd.melt(pi,id_vars=["start","end"],value_name="pi",var_name="sample")



covplot = ggplot(covpcm, aes(x="sample",y="covpc"))+ geom_violin() +coord_flip()
covplot.save((outfile+"_covpc.png"),width=150,height=150,units="mm")

dpplot = ggplot(meandepthm, aes(x="sample",y="depth"))+ geom_violin() +coord_flip()
dpplot.save((outfile+"_meandp.png"),width=150,height=150,units="mm")

dpplot = ggplot(pim, aes(x="sample",y="pi"))+ geom_violin() +coord_flip()
dpplot.save((outfile+"_pi.png"),width=150,height=150,units="mm")
