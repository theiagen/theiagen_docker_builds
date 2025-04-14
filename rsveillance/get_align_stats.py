#!/usr/bin/env python
import argparse

# Set up argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--sample', '-s', required=True, help="sample name")
parser.add_argument('--target', '-t', required=True, help="sarget name")
parser.add_argument('--subfactor', '-F', type=int, required=True, help="depth correction factor (subsampling proportion)")
parser.add_argument('--flagstats', '-i', required=True, help="path to the flagstats file")
parser.add_argument('--dhist', '-d', required=True, help="path to depth histogram file")
parser.add_argument('--output', '-o', required=True, help="outfile prefix")
parser.add_argument('--mindepth', '-m', type=int, default=10, help="minimum depth")

args = parser.parse_args()

def process_stats(sample, target, subfact, flagstats_file, dhist_file, output_file, mindepth):
    #get reads aligned
    reads = -1
    aligned = -1
    paired = -1
    with open(flagstats_file, "r") as f:
        for l in f:
            l = l.strip().split('\t')
            passreads = l[0]
            failreads = l[1]
            stat = l[2]
            if stat == 'total (QC-passed reads + QC-failed reads)': 
                reads = int(passreads)
            elif stat == "mapped": 
                aligned = int(passreads)
            elif stat == "properly paired": 
                paired = int(passreads)

    #get coverage / depth
    goodcov = 0
    cov = 0
    gsize = 0
    dtotal = 0
    with open(dhist_file, "r") as f:
        for l in f:
            l = l.split("\t")
            depth = int(l[1])
            count = int(l[2])
            cdepth = int(l[3])
            dtotal += cdepth * count
            if cdepth >= mindepth:
                gsize += count
                cov += count
                goodcov += count
            elif cdepth > 0:
                gsize += count
                cov += count
            elif cdepth == 0:
                gsize += count

    meandepth = round(dtotal / gsize)

    # Write stats to the output file
    with open(output_file, "w") as f:
        print("\t".join(map(str, [sample, target, subfact,
                                  reads, aligned, paired,
                                  meandepth, goodcov, cov, gsize])), file=f)

output_file = args.output + "_alignstats.txt"

process_stats(args.sample, args.target, args.subfactor, args.flagstats, args.dhist, output_file, args.mindepth)
