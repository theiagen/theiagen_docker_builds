#!/usr/bin/env python

import numpy as np
import sys
from Bio import SeqIO, Seq

import pandas as pd

import re
import sys
import os

import argparse


mingoodcov = 0.8
minratio = 0.95

parser = argparse.ArgumentParser()
parser.add_argument('--alignstats', '-i', help='align stats file')
parser.add_argument('--coverage', '-c', help='min coverage for secondary call')
parser.add_argument('--ratio', '-r', help='min ratio for secondary call')
parser.add_argument('--out', '-o', help='outfile')

args = parser.parse_args()
alignstatsf = args.alignstats
mingoodcov = float(args.coverage)
minratio = float(args.ratio)
outfile = args.out


#read in alignment stats
alignstatsnames = ["sample","mashcall","subfactor","reads","aligned","paired","meandepth","goodcov","cov","gsize"]
alignstats = pd.read_csv(alignstatsf,sep="\t",names=alignstatsnames)

#add coverage percentage to alignstats
alignstats['covpc'] = alignstats['cov'] / alignstats['gsize']
alignstats['goodcovpc'] = alignstats['goodcov'] / alignstats['gsize']


#filter low pc coverage aligns
failedcov = len(alignstats[alignstats['goodcovpc'] < mingoodcov])
print("removing {}/{} calls due to coverage below {}".format(failedcov,len(alignstats),mingoodcov),
      file=sys.stderr)
alignstats = alignstats[alignstats['goodcovpc'] > mingoodcov]


#identify single samples (align to only RSV A or B)
singles = alignstats.loc[alignstats['sample'].duplicated(keep=False)==False]

#identify duplicated samples (align to both RSV A and B)
dupes = alignstats.loc[alignstats['sample'].duplicated(keep=False)]

#if any duplicates are found, investigate to see which are true mixes
if np.shape(dupes)[0] >=1:
    dupescf = dupes.pivot(index="sample",
                            values=["meandepth","covpc","goodcovpc"],
                            columns="mashcall",)
    dupescf["goodcovratio"] = dupescf["goodcovpc"]["RSVA"]/dupescf["goodcovpc"]["RSVB"]
    dupescf["meandepthratio"] = dupescf["meandepth"]["RSVA"]/dupescf["meandepth"]["RSVB"]


    #identify putative coinfections (second align is within minratio of first)
    coinfections = (dupescf["goodcovratio"] > minratio) & (dupescf["goodcovratio"] < (1/minratio))
    rsva = (dupescf["goodcovratio"] > (1/minratio))
    rsvb = (dupescf["goodcovratio"] < minratio)

    #assemble final calls table
    calls = pd.concat([pd.DataFrame({"sample":dupescf.index[rsvb].tolist(),
                            "call":["RSVB"]*sum(rsvb)}),
                        pd.DataFrame({"sample":dupescf.index[rsva].tolist(),
                            "call":["RSVA"]*sum(rsva)}),
                        pd.DataFrame({"sample":dupescf.index[coinfections].tolist(),
                            "call":["RSVA/B"]*sum(coinfections)}),
                        singles[['sample',"mashcall"]]
                        ], axis=0)

    #combine mash / final calls into one table
    combinations = pd.concat([pd.DataFrame({"sample":dupescf.index[rsvb].tolist(),
                            "mashcall":["RSVB"]*sum(rsvb),
                            "call":["RSVB"]*sum(rsvb)}),
                            pd.DataFrame({"sample":dupescf.index[rsva].tolist(),
                                            "mashcall":["RSVA"]*sum(rsva),
                                            "call":["RSVA"]*sum(rsva)}),
                            pd.DataFrame({"sample":dupescf.index[coinfections].tolist(),
                                            "mashcall":["RSVB"]*sum(coinfections),
                                            "call":["RSVA/B"]*sum(coinfections)}),
                            pd.DataFrame({"sample":dupescf.index[coinfections].tolist(),
                                            "mashcall":["RSVA"]*sum(coinfections),
                                            "call":["RSVA/B"]*sum(coinfections)}),
                            pd.DataFrame({"sample":singles['sample'].to_list(),
                                            "mashcall":singles['mashcall'].to_list(),
                                            "call":singles['mashcall'].to_list()}),
                            ], axis=0)
    
else:
    #use single A/B calls as final calls table and print to file
    calls = singles[['sample',"mashcall"]]
    combinations = pd.DataFrame({"sample":singles['sample'].to_list(),
                "mashcall":singles['mashcall'].to_list(),
                "call":singles['mashcall'].to_list()})

#merge calls with alignstats
combinations = pd.merge(combinations,alignstats[alignstatsnames],on=['sample','mashcall'])

calls.to_csv("{}_calls.txt".format(outfile),sep="\t",index=False)
combinations.to_csv("{}_alignstats.txt".format(outfile),sep="\t",index=False,header=False)
