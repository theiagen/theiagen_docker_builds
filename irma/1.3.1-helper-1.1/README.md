# IRMA Image and Helper Script

This helper script is used to perform post analysis reformatting of FASTA files and QC metrics gathering. The script will ascertain the subtype of HA and NA segments, create padded and concatenated consensus files, and create a MIRA-like QC summary file.

The Docker image associated with this script utilizes the Staphb IRMA image as base and adds pertinent Python packages to perform the FASTA manipulation, writing, and QC table creation.

- staphb/irma:1.3.1
- biopython
- pandas

### Running the Script
A bare bones testing set has been included under `test_data/flu_test`. For this particular set of data a test can be run as follows. This usage mimics how the script it utilized within the IRMA WDL task. Info regarding missing variant and coverage tables is expected behavior for this dataset as only A_HA_H3 variant tables were included. The resulting QC table will have only the HA_H3 segment populated and the rest should be "N/A".
```
python irma_helper.py -d test_data/flu_test/ -t A -s flu_test -v
```
A Flu B set is included under `test_data/fluB_test`. This is IRMA 1.3.1 output from simulated paired-end reads (generated from the IRMA FLU module's B references with a few SNPs and a ~20% minor variant at B_HA position 901), so it does not represent a real sample. Only the B_HA tables and references were kept, so the resulting QC table will have coverage and variant metrics for HA only and "N/A" for the rest; read counts are populated for all segments. Segment 1 should be labeled PB1 and segment 2 PB2.
```
python irma_helper.py -d test_data/fluB_test/ -t B -s fluB_test -v
```
Note that the helper renames the files in `amended_consensus/` and creates `padded_assemblies/` next to the input directory, so run it against a copy of the test data.
```bash
usage: irma_helper.py [-h] -d INPUT_DIR -t IRMA_TYPE -s SAMPLENAME [--log LOG] [-v]

A helper script to parse IRMA output files and create various consensus FASTA files

options:
  -h, --help            show this help message and exit
  -d INPUT_DIR, --input_dir INPUT_DIR
                        Output directory of IRMA used as input for helper script (default: None)
  -t IRMA_TYPE, --irma_type IRMA_TYPE
                        IRMA Flu Type; A or B (default: None)
  -s SAMPLENAME, --samplename SAMPLENAME
                        Samplename used to build filepaths (default: None)
  --log LOG             Log file (default: irma_helper.log)
  -v, --verbose         Verbose output (default: False)
```

## Changelog

- 1.3.1-helper-1.1 fixes the following bugs in the 1.3.1-helper script:
  - The Flu B segment order swapped PB1 and PB2 (segment 1 was labeled PB2 and segment 2 was labeled PB1), so those two segments were mislabeled in the consensus FASTA filenames and headers (including the padded copies) for Type B samples. Segment 1 is now PB1 and segment 2 is PB2.
  - The unpadded multi-FASTA (`amended_consensus/<samplename>.irma.consensus.fasta`) contained padded sequences (periods replaced by Ns and dashes removed) because the padding step modified the same sequence record that was later written to the unpadded file.
  - Type A samples without an HA or NA subtype used "No subtype" as a filename suffix, so the coverage, reference, and variant files for those segments were never found and the QC summary reported `HA_No subtype` with N/A values.
  - Total Reads and Pass QC were reported as N/A for any segment missing from `READ_COUNTS.tsv`, even though these are sample-level counts. A missing `READ_COUNTS.tsv` file also crashed the script instead of reporting N/A.
  - A segment with a coverage file but no reference file crashed the script when calculating `% Reference Covered`; this is now reported as N/A.
  - The concatenated variant tables (`<samplename>_irma_all_{variants,insertions,deletions}.tsv`) rounded frequencies and rewrote `TRUE`/`FALSE` as `True`/`False`. Rows are now copied exactly as IRMA wrote them.
