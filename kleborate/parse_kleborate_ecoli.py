#!/usr/bin/env python3
import csv
import sys
import argparse

def parse_kleborate_ecoli_output(input_file):
    """Parse Kleborate TSV output and write results to individual files."""
    
    def clean_value(value):
        """Convert '-' to empty string, otherwise return the value. Keep list in case need to expand"""
        return '' if value in ['-'] else value
    
    with open(input_file, 'r') as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        tsv_data = list(tsv_reader)
        tsv_dict = dict(zip(tsv_data[0], tsv_data[1]))
        
        # Write species
        with open("SPECIES", 'wt') as Species:
            kleb_species = clean_value(tsv_dict.get('species', ''))
            Species.write(kleb_species)
            
        # Write LEE MLST sequence type
        with open("LEE_ST", 'wt') as LEE_MLST_Sequence_Type:
            mlst = clean_value(tsv_dict.get('LEE_ST', ''))
            LEE_MLST_Sequence_Type.write(mlst)

        # Write LEE lineage
        with open("LEE_LINEAGE", 'wt') as LEE_Lineage:
            lineage = clean_value(tsv_dict.get('LEE_lineage', ''))
            LEE_Lineage.write(lineage)
            
        # Write Pathotype
        with open("PATHOTYPE", 'wt') as Pathotype:
            pathotype = clean_value(tsv_dict.get('Pathotype', ''))
            Pathotype.write(pathotype)


def main():
    parser = argparse.ArgumentParser(
        description='Parse Kleborate TSV output and extract key information'
    )
    parser.add_argument(
        'input_file',
        help='Path to the Kleborate TSV output file'
    )
    
    args = parser.parse_args()
    
    try:
        parse_kleborate_ecoli_output(args.input_file)
        print(f"Successfully parsed {args.input_file}")
    except FileNotFoundError:
        print(f"Error: File '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()