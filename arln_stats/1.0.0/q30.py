#!/usr/bin/env python3

import os
import sys
import gzip
import argparse

# Disable cache usage in Python so __pycache__ isn't formed.
sys.dont_write_bytecode = True

# Function to get the script version
def get_version():
    return "2.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version', version=get_version())  # Add an argument to display the version
    parser.add_argument('-i', '--input_reads', dest='input', required=True, help='FASTQs', nargs='+')
    return parser.parse_args()

def qual_stat(qstr):
    q20 = 0
    q30 = 0
    for q in qstr:
        qual = ord(q) - 33  # Remove chr() to directly use the character
        if qual >= 30:
            q30 += 1
            q20 += 1
        elif qual >= 20:
            q20 += 1
    return q20, q30

def stat(filename, read_num):
    total_read_count = 0
    total_base_count = 0
    q20_count = 0
    q30_count = 0

    # Check if the input file is gzipped
    open_func = gzip.open if filename.endswith('.gz') else open

    with open_func(filename, 'rt') as f:
        while True:
            header = f.readline()
            if not header:  # End of file
                break
            sequence = f.readline().strip()
            plus_line = f.readline()  # Skip the '+' line
            quality_str = f.readline().strip()

            total_read_count += 1
            total_base_count += len(sequence)
            q20, q30 = qual_stat(quality_str)
            q20_count += q20
            q30_count += q30
    
    print(read_num, 100 * (float(q30_count) / float(total_base_count))) 

def main():
    args = parseArgs()

    for i, filename in enumerate(args.input):
        read_num = f"read{i+1}"
        stat(filename, read_num)

if __name__ == "__main__":
    main()