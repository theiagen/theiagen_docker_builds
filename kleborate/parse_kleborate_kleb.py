#!/usr/bin/env python3
import csv
import sys
import argparse

def parse_kleborate_output(input_file):
    """Parse Kleborate TSV output and write results to individual files."""
    
    with open(input_file, 'r') as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        tsv_data = list(tsv_reader)
        tsv_dict = dict(zip(tsv_data[0], tsv_data[1]))
        
        # Write species
        with open("SPECIES", 'wt') as Species:
            kleb_species = tsv_dict.get('species', '')
            Species.write(kleb_species)
        
        # Write MLST sequence type
        with open("MLST_SEQUENCE_TYPE", 'wt') as MLST_Sequence_Type:
            mlst = tsv_dict.get('ST', '')
            MLST_Sequence_Type.write(mlst)
        
        # Write virulence score
        with open("VIRULENCE_SCORE", 'wt') as Virulence_Score:
            virulence_level = tsv_dict.get('virulence_score', '')
            Virulence_Score.write(virulence_level)
        
        # Write resistance score
        with open("RESISTANCE_SCORE", 'wt') as Resistance_Score:
            resistance_level = tsv_dict.get('resistance_score', '')
            Resistance_Score.write(resistance_level)
        
        # Write number of resistance genes
        with open("NUM_RESISTANCE_GENES", 'wt') as Num_Resistance_Genes:
            resistance_genes_count = tsv_dict.get('num_resistance_genes', '')
            Num_Resistance_Genes.write(resistance_genes_count)
        
        # Write BLA resistance genes
        with open("BLA_RESISTANCE_GENES", 'wt') as BLA_Resistance_Genes:
            bla_res_genes_list = ['Bla_acquired', 'Bla_inhR_acquired', 'Bla_ESBL_acquired', 
                                  'Bla_ESBL_inhR_acquired', 'Bla_Carb_acquired']
            bla_res_genes = []
            for i in bla_res_genes_list:
                value = tsv_dict.get(i, '-')
                if value != '-':
                    bla_res_genes.append(value)
            bla_res_genes_string = ';'.join(bla_res_genes)
            BLA_Resistance_Genes.write(bla_res_genes_string)
        
        # Write ESBL resistance genes
        with open("ESBL_RESISTANCE_GENES", 'wt') as ESBL_Resistance_Genes:
            esbl_res_genes_list = ['Bla_ESBL_acquired', 'Bla_ESBL_inhR_acquired']
            esbl_res_genes = []
            for i in esbl_res_genes_list:
                value = tsv_dict.get(i, '-')
                if value != '-':
                    esbl_res_genes.append(value)
            esbl_res_genes_string = ';'.join(esbl_res_genes)
            ESBL_Resistance_Genes.write(esbl_res_genes_string)
        
        # Write key resistance genes
        with open("KEY_RESISTANCE_GENES", 'wt') as Key_Resistance_Genes:
            key_res_genes_list = ['Col_acquired', 'Fcyn_acquired', 'Flq_acquired', 'Rif_acquired', 
                                  'Bla_acquired', 'Bla_inhR_acquired', 'Bla_ESBL_acquired', 
                                  'Bla_ESBL_inhR_acquired', 'Bla_Carb_acquired']
            key_res_genes = []
            for i in key_res_genes_list:
                value = tsv_dict.get(i, '-')
                if value != '-':
                    key_res_genes.append(value)
            key_res_genes_string = ';'.join(key_res_genes)
            Key_Resistance_Genes.write(key_res_genes_string)
        
        # Write genomic resistance mutations
        with open("GENOMIC_RESISTANCE_MUTATIONS", 'wt') as Resistance_Mutations:
            res_mutations_list = ['Bla_chr', 'SHV_mutations', 'Omp_mutations', 
                                  'Col_mutations', 'Flq_mutations']
            res_mutations = []
            for i in res_mutations_list:
                value = tsv_dict.get(i, '-')
                if value != '-':
                    res_mutations.append(value)
            res_mutations_string = ';'.join(res_mutations)
            Resistance_Mutations.write(res_mutations_string)
        
        # Write K type
        with open("K_TYPE", 'wt') as Ktype:
            ktype = tsv_dict.get('K_type', '')
            Ktype.write(ktype)
        
        # Write K locus
        with open("K_LOCUS", 'wt') as Klocus:
            klocus = tsv_dict.get('K_locus', '')
            Klocus.write(klocus)
        
        # Write O type
        with open("O_TYPE", 'wt') as Otype:
            otype = tsv_dict.get('O_type', '')
            Otype.write(otype)
        
        # Write O locus
        with open("O_LOCUS", 'wt') as Olocus:
            olocus = tsv_dict.get('O_locus', '')
            Olocus.write(olocus)
        
        # Write K locus confidence
        with open("K_LOCUS_CONFIDENCE", 'wt') as K_locus_confidence:
            k_locus_confidence = tsv_dict.get('K_locus_confidence', '')
            K_locus_confidence.write(k_locus_confidence)
        
        # Write O locus confidence
        with open("O_LOCUS_CONFIDENCE", 'wt') as O_locus_confidence:
            o_locus_confidence = tsv_dict.get('O_locus_confidence', '')
            O_locus_confidence.write(o_locus_confidence)

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
        parse_kleborate_output(args.input_file)
        print(f"Successfully parsed {args.input_file}")
    except FileNotFoundError:
        print(f"Error: File '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()