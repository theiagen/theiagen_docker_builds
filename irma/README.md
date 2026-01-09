# IRMA Image and Helper Script

This helper script is used to perform post analysis reformatting of FASTA files and QC metrics gathering. The script will ascertain the subtype of HA and NA segments, create padded and concatenated consensus files, and create a MIRA-like QC summary file. 

The Docker image associated with this script utilizes the Staphb IRMA image as base and adds pertinent Python packages to perform the FASTA manipulation, writing, and QC table creation. 

- staphb/irma:1.3.1
- biopython
- pandas

### Running the Script
A bare bones testing set has been included under `test_data/flu_test`. The script is designed to be run within `test_data`. For this particular set of data a test can be run as follows. This usage mimics how the script it utilized within the IRMA WDL task. 
```
../irma_helper.py -d flu_test -t A -s flu_test -v
```
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