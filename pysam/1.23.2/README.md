# Pysam container

Main tool: [Pysam](https://github.com/pysam-developers/pysam) | [samtools](https://github.com/samtools/samtools)

Code repository: [theiagen_docker_builds](https://github.com/theiagen/theiagen_docker_builds)

Additional tools:
- biopython 1.87
- numpy 2.2.6
- pysam 0.24.0
- samtools 1.23 
- contaminant_check.py NA
- gene_coverage.py NA
- variant_annotation.py NA

## Changes
Add `variant_annotation.py` to translate query genes with detected variants and
annotate each variant's protein-level consequence (missense/synonymous/nonsense
substitutions, in-frame insertions/deletions and frameshifts).  It is run
separately from `gene_coverage.py` and takes a VCF (either a raw VCF or a
`GENE_VARIANTS.vcf` produced by `gene_coverage.py`) together with a reference
`--reference_gbff` (or `--reference_gff` and `--reference_fa`), writing the
report to `VARIANT_ANNOTATIONS.txt`.

## Example Usage

### Contaminant Check
```
python3 /usr/bin/contaminant_check.py -h
usage: contaminant_check.py [-h] --expected_sequences EXPECTED_SEQUENCES
                            [EXPECTED_SEQUENCES ...]
                            [--minimum_expected_sequences MINIMUM_EXPECTED_SEQUENCES]
                            [--maximum_unexpected_sequences MAXIMUM_UNEXPECTED_SEQUENCES]
                            --coverage_by_sequence_json COVERAGE_BY_SEQUENCE_JSON
                            --depth_by_sequence_json DEPTH_BY_SEQUENCE_JSON
                            --reads_by_sequence_json READS_BY_SEQUENCE_JSON
                            --contaminant_fasta CONTAMINANT_FASTA --minimum_percent_coverage
                            MINIMUM_PERCENT_COVERAGE --minimum_depth MINIMUM_DEPTH
                            --minimum_reads_mapped MINIMUM_READS_MAPPED
```

### Gene Coverage
```
python3 /usr/bin/gene_coverage.py -h  
usage: gene_coverage.py [-h] --bam BAM [--bedfile BEDFILE] [--reference_gbff REFERENCE_GBFF]
                        --query_genes QUERY_GENES [QUERY_GENES ...]
                        [--feature_type FEATURE_TYPE]
                        [--feature_qualifier FEATURE_QUALIFIER] [--exact_match]
                        [--min_depth MIN_DEPTH]
```

### Variant Annotation
```
python3 /usr/bin/variant_annotation.py -h
usage: variant_annotation.py [-h] --vcf VCF --reference_gbff REFERENCE_GBFF
                             --query_genes QUERY_GENES [QUERY_GENES ...]
                             [--feature_type FEATURE_TYPE]
                             [--feature_qualifier FEATURE_QUALIFIER] [--exact_match]
                             [--transl_table TRANSL_TABLE] [--output OUTPUT]
```

Example report:
```
lanosterol.14-alpha.demethylase: lanosterol 14-alpha demethylase (missense_variant c.395A>T p.Tyr132Phe; A:20 T:0)
```

