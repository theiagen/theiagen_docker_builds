#!/usr/bin/env python

import re
import sys
import argparse 

from statistics import mean

# vcf='pileup_all_vars.vcf.gz'
# ref="reference.fasta"
# #gff="reference.gff"
# rangesbed="TBscheme.insert.bed"
# depthfile="all_depth.txt.gz"
# outfile="TB001_amplicon"

winsize = 1000

parser = argparse.ArgumentParser()
parser.add_argument('--window', '-w', help='depth window size')
parser.add_argument('--depth', '-d', help='depth file (samtools / txt)')
parser.add_argument('--out', '-o', help='outfile prefix')
parser.add_argument('--sample', '-s', help='outfile prefix')
parser.add_argument('--depthfactor', '-F', help='depth correction factor (subsampling proportion)')

args = parser.parse_args()
winsize = int(args.window)
depthfile = args.depth
depthfactor = int(args.depthfactor)
sample = args.sample
outfile = args.out


def printblock(posn,dep,outfile):
    if sample and depthfactor:
        print("{}\t{}\t{}\t{}".format(sample,posn,round(dep),round(dep*depthfactor)), file=outfile)
    elif depthfactor:
        print("{}\t{}\t{}".format(posn,round(dep),round(dep*depthfactor)), file=outfile)
    elif sample:
        print("{}\t{}\t{}".format(sample,posn,round(dep)), file=outfile)
    else:
        print("{}\t{}".format(posn,round(dep)), file=winsfile)

def printhist(dep,count,outfile):
    if sample and depthfactor:
        print("{}\t{}\t{}\t{}".format(sample,round(dep),count,round(dep*depthfactor)), file=outfile)
    elif depthfactor:
        print("{}\t{}\t{}".format(round(dep),count,round(dep*depthfactor)), file=outfile)
    elif sample:
        print("{}\t{}\t{}".format(sample,round(dep),count), file=outfile)
    else:
        print("{}\t{}".format(round(dep),count), file=winsfile)

dhist = dict()
windepths=list()

winsfile = open(outfile+"_depthwins.txt",'w')
with open(depthfile) as my_file:
    for line in my_file:
        line = line.split()
        c = line[0]
        if c == "#CHROM": continue
        s = int(line[1])
        d = int(line[2])
        if d not in dhist:
            dhist[d]=0
        dhist[d] +=1
        windepths.append(d)
        #print("{}  {}  {}".format(s,winsize,s%winsize),file=sys.stdout)
        if s % winsize == 0:
            printblock(s,mean(windepths),winsfile)
            windepths=list()

if s % winsize > 0:
    printblock(s,mean(windepths),winsfile)
winsfile.close()


histfile = open(outfile+"_depthhist.txt",'w')
for d in sorted(dhist.keys()):
    #if d not in dhist:
    #    dhist[d]=0
    printhist(d,dhist[d],histfile)
histfile.close()
