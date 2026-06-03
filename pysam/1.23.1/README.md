# Pysam container

Main tool: [Pysam](https://github.com/pysam-developers/pysam) | [samtools](https://github.com/samtools/samtools)

Code repository: [theiagen_docker_builds](https://github.com/theiagen_docker_builds/)

Additional tools:
- biopython 1.87
- numpy 2.2.6
- pysam 0.24.0
- samtools 1.23 
- contaminant_check.py NA
- gene_coverage.py NA

## Changes
Add `gene_coverage.py` and `biopython`

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
