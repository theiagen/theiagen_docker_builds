#!/usr/bin/env python

import allel
import numpy as np

import re
import sys
import argparse 
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--bed', '-b', help='primer / amplicon positions (bed)')
parser.add_argument('--gff', '-g', help='gene positions (gff)')
parser.add_argument('--depth', '-d', help='depth file (samtools / txt)')
parser.add_argument('--sample', '-s', help='sample name')
parser.add_argument('--factor', '-F', help='subsampling factor')
parser.add_argument('--out', '-o', help='outfile prefix')

args = parser.parse_args()
bedfile = args.bed
gfffile = args.gff
depthfile = args.depth
sample = args.sample
outfile = args.out


depth = pd.read_csv(depthfile, sep='\t')

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

print(ranges)
nranges=len(ranges)

meandepth=pd.DataFrame({"sample":[sample]*nranges,
              "name":[n for s,e,n in ranges],
              "start":[s for s,e,n in ranges],
              "end":[e for s,e,n in ranges],
              "depth":[-1]*nranges
              })

i=-1
for st,en,name in ranges:
    i=i+1
    indepth = np.logical_and((depth['POS'] >= st),(depth['POS'] <= en))

    if sum(indepth)>0:
        locdepth = round(depth[indepth].iloc[:,2].mean().tolist())
    else:
        locdepth = -1
    meandepth.iloc[i,:] = [sample,name,st,en,locdepth]


meandepth.to_csv(outfile,sep="\t",index=False,)