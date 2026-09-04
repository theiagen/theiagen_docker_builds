# Pysam container

Main tool: [Pysam](https://github.com/pysam-developers/pysam) | [samtools](https://github.com/samtools/samtools)

Code repository: [theiagen_docker_builds](https://github.com/theiagen/theiagen_docker_builds)

Additional tools:
- pysam 0.24.0
- samtools 1.23
- contaminant_check.py NA

## Changes
Migrate `gene_coverage.py` and `variant_annotation.py` out of this container into
the standalone [`theiagene`](https://github.com/theiagen/theiagene) package.
`biopython` was only required by those two tools and has been dropped from this
image.

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
