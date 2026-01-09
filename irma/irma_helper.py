#!/usr/bin/env python3
import sys
import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from typing import Dict
"""
A helper script to facilitate QC summary creation and FASTA concatenation.
"""

# flu segments from largest to smallest
# Thank you Molly H. for this code block!
# declare associative arrays for segment numbers
FLU_SEGMENTS: Dict[str, Dict[str, str]] = {
"A": {"1": "PB2", "2": "PB1", "3": "PA", "4": "HA", "5": "NP", "6": "NA", "7": "MP", "8": "NS"},
"B": {"2": "PB1", "1": "PB2", "3": "PA", "4": "HA", "5": "NP", "6": "NA", "7": "MP", "8": "NS"},
}

def setup_logger(log_file: str, verbose: bool) -> logging.Logger:
  log_level = logging.DEBUG if verbose else logging.INFO
  
  logger = logging.getLogger(__name__)
  logger.setLevel(log_level)
  
  file_handler = logging.FileHandler(log_file)
  stream_handler = logging.StreamHandler(sys.stderr)

  formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")
  file_handler.setFormatter(formatter)
  stream_handler.setFormatter(formatter)
  
  logger.addHandler(file_handler)
  logger.addHandler(stream_handler)
  
  return logger

def subtype_selection(
                      input_dir: Path,
                      samplename: str,
                      irma_type: str,
                      logger: logging.Logger
):
  subtype_dict = {"HA" : "", "NA" : ""}
  subtype_notes = ""

  # Defensively assign flu segments 
  if irma_type not in FLU_SEGMENTS:
    raise ValueError("IRMA Type not of type A or B")
  flu_segments = FLU_SEGMENTS[irma_type]

  if irma_type == "A":
    logger.info(f"Type reported as: {irma_type}")
    
    # Check for any subtypes within Type A segments, relying on structure of A_{segment}_{subtype}.fasta
    # Search is filtered to HA and NA segments only in the top layer of directory.
    for file in os.scandir(input_dir):
      if file.is_file():
        # Splitting these two into different paths preserves the order of the subtype 'H*N*'
        if file.name.startswith("A_HA_H") and file.name.endswith(".fasta"):
          logger.info(f"Aquiring HA subtype from {file}.")
          subtype_dict["HA"] = os.path.splitext(file.name)[0].split("_")[-1]
        elif file.name.startswith("A_NA_N") and file.name.endswith(".fasta"):
          logger.info(f"Aquiring NA subtype from {file}.")
          subtype_dict["NA"] = os.path.splitext(file.name)[0].split("_")[-1]
      else:
        logger.info(f"File {file.name} does not match segment criteria or does not contain subtype information.")

    # Build subtype in correct order
    subtype = f"{subtype_dict['HA']}{subtype_dict['NA']}"
    if subtype != "":
      logger.info(f"Subtype reported as {subtype_dict}")
    else:
      logger.info(f"No subtype found within {samplename} for type {irma_type}.")
      subtype_dict = {"HA": "No subtype", "NA": "No subtype"}
      subtype = "No subtype predicted by IRMA"

  # Populate subtype and subtype notes with set values for Type B
  elif irma_type == "B":
    logger.info(f"Type reported as: {irma_type}")
    subtype = "No subtype predicted by IRMA"
    subtype_notes = "IRMA does not differentiate Victoria and Yamagata Flu B lineages. See abricate_flu_subtype output column"
    logger.info(f"IRMA does not differentiate Victoria and Yamagata Flu B lineages. See abricate_flu_subtype output column")
  
  # Write subtype and subtype notes to file for WDL interpretation
  try:
    logger.debug("Writing IRMA_SUBTYPE to file.")
    with open("IRMA_SUBTYPE.txt", "w") as irma_subtype:
      irma_subtype.writelines(subtype)
  except IOError as e:
    logger.error(f"Error writing IRMA_SUBTYPE.txt file: {e}")

  try:
    logger.info("Writing IRMA_SUBTYPE_NOTES to file.")
    with open("IRMA_SUBTYPE_NOTES.txt", "w") as subtype_notes_file:
      subtype_notes_file.writelines(subtype_notes)
  except IOError as e:
    logger.error(f"Error writing IRMA_SUBTYPE_NOTES.txt file: {e}")

  return flu_segments, subtype_dict

def consensus_creation(
                        input_dir: Path,
                        samplename: str,
                        flu_segments: Dict[str, str],
                        logger: logging.Logger
):
  
  consensus_array = []
  padded_consensus_array = []
  concatenated_seq = ""

  # Build consensus files using flu_segments as a guide
  for segment_idx in flu_segments:
    segment = flu_segments[segment_idx]
    logger.debug(f"Parsing {segment} in output files")
    segment_file = input_dir / f"amended_consensus/{samplename}_{segment_idx}.fa"
    # Check if the segment file is present, it is possible for it not to be
    if segment_file.exists(): 
      logger.debug(f"Segment file selected: {segment_file}")
      # Create a SeqIO record and rename FASTA
      record = list(SeqIO.parse(segment_file, "fasta"))[0]
      record.id = record.id.replace(segment_idx, segment)
      record.description = ""

      concatenated_seq += str(record.seq)
      consensus_array.append(record)

      # Write FASTA
      try:
        SeqIO.write(record, segment_file, "fasta-2line")
        segment_file.rename(input_dir / f"amended_consensus/{samplename}_{segment_idx}.fa")
        logger.debug(f"Amended consensus FASTA for {samplename}_{segment} written and renamed")
      except IOError as e:
        logger.error(f"Error writing segment fasta: {e}")

      # Write padded FASTA
      try:
        # Replace . for N and remove -
        record.seq = Seq(str(record.seq).replace('.', 'N').replace('-', ''))
        SeqIO.write(record, f"padded_assemblies/{samplename}_{segment}.pad.fasta", "fasta-2line")
        logger.debug(f"Padded fasta written for: {segment_file}")
      except IOError as e:
        logger.error(f"Error writing padded segment fasta: {e}")

      padded_consensus_array.append(record)
    else:
      logger.warning(f"No file found for segment: {segment}")
      continue
  
  # Write array of consensus FASTAs
  try:
    SeqIO.write(consensus_array, input_dir / f"amended_consensus/{samplename}.irma.consensus.fasta", "fasta-2line")
    SeqIO.write(padded_consensus_array, f"padded_assemblies/{samplename}.irma.consensus.pad.fasta", "fasta-2line")
    SeqIO.write(SeqRecord(Seq(concatenated_seq), f"{samplename}_irma_concatenated", description=""), input_dir / f"amended_consensus/{samplename}.irma.consensus.concatenated.fasta", "fasta-2line")
    SeqIO.write(SeqRecord(Seq(str(concatenated_seq).replace('.', 'N').replace('-', '')), f"{samplename}_irma_concatenated_padded", description=""), f"padded_assemblies/{samplename}.irma.consensus.concatenated.pad.fasta", "fasta-2line")
  except IOError as e:
    logger.error(f"Error writing FASTA files: {e}")

def read_counts(
                input_dir: Path,
                samplename: str, 
                logger: logging.Logger, 
                segment: str
):
  # Read in READ_COUNTS.tsv and pull read metrics
  try:
    read_counts = pd.read_csv( input_dir / f"tables/READ_COUNTS.tsv", sep='\t', index_col=0)
    logger.debug(f"READ_COUNTS.tsv file found for {samplename}.")

    # Obtain mapped reads per segment
    if not read_counts.empty:
      logger.debug("READ_COUNTS.tsv present, parsing mapped reads")
      reads_mapped = read_counts.loc[read_counts.index.str.contains(segment), 'Reads'].values[0]
      logger.info(f"Mapped reads found to be: {reads_mapped}")

    total_reads = read_counts.loc["1-initial", "Reads"]
    pass_qc_reads = read_counts.loc["2-passQC", "Reads"]
    logger.info(f"Total Reads found to be: {total_reads} \n Passing Reads found to be: {pass_qc_reads}")

  except FileNotFoundError:
    logger.warning(f"WARNING: READ_COUNTS.tsv file not found for {samplename}. Cannot extract read counts for QC summary.")
  return total_reads, pass_qc_reads, reads_mapped

def reference_length(
                    samplename: str, 
                    logger: logging.Logger, 
                    segment: str, 
                    reference_path: Path
):
  # Load reference seq to obtain length of seq
  try:
    ref_record = list(SeqIO.parse(reference_path, "fasta"))[0]
    segment_ref_len=str(len(ref_record))
    logger.debug(f"Segment reference found with length {segment_ref_len} for {samplename}:{segment}")
  except FileNotFoundError:
    logger.warning(f"WARNING: No reference file found for segment {segment} for {samplename}")
    segment_ref_len = "N/A"
  return segment_ref_len

def coverage_parsing(samplename: str,
                     logger: logging.Logger, 
                     segment: str, 
                     ref_length: int, 
                     coverage_path: Path
):
  # Load coverage file for parsing
  try:
    coverages = pd.read_csv(coverage_path, sep='\t')
    logger.debug(f"Coverage file found for {samplename}:{segment}")
    mapped_bases = 0
    # Mapped bases
    for base in coverages['Consensus']:
      if base not in ['-','N', 'a', 't', 'c', 'g']:
        mapped_bases += 1
    # Calculate percentage rounded to second decimal place
    if mapped_bases > 0 and int(ref_length) > 0:
      segment_pct_ref_cov = round((mapped_bases / int(ref_length) * 100), 2)
    else:
      segment_pct_ref_cov = 0
    # Obtain mean and median coverage from coverage file
    mean_cov = round(coverages['Coverage Depth'].mean(), 2)
    median_cov = round(coverages['Coverage Depth'].median(), 2)
    logger.info(f"Mean Coverage and Median Coverage calculated as {mean_cov} and {median_cov} respectively.")
  except FileNotFoundError as e:
    logger.warning(f"Error coverage file not found for {samplename}:{segment}, setting mean_cov and median_cov to 'N/A'")
    mean_cov = "N/A"
    median_cov = "N/A"
    segment_pct_ref_cov = None
  return mean_cov, median_cov, segment_pct_ref_cov

def variant_parsing(logger: logging.Logger,
                    variant_location: Path
):
  var_types = ["variants", "insertions", "deletions"]
  variant_threshold = 0.05

  snv_files = pd.DataFrame()
  deletion_files = pd.DataFrame()
  insertion_files = pd.DataFrame()

  # Set default variant counts
  variant_count = ""
  insertion_count = ""
  deletion_count = ""

  # Using var_types reference dictionary obtain variant frequency for each type per segment
  # and append to concatenated variant files
  for variant in var_types:
    logger.debug(f"Checking type {variant}")
    logger.debug(f"Variant location: {variant_location}")
    variant_file_path = Path(f"{variant_location}-{variant}.txt")

    logger.debug(f"Variant file selected: {variant_file_path}")
    try:
      variant_df = pd.read_csv(variant_file_path, sep='\t')
    except FileNotFoundError:
      logger.error(f"Variant not found for variant: {variant}")
      continue
    logger.debug("Variant file loaded as dataframe")

    try: 
      # Perform frequency count on variants file
      if variant == "variants":
        logger.debug(f"Variant file of type '{variant}' selected")
        variant_count = str((variant_df["Minority_Frequency"] >= variant_threshold).sum())
        logger.debug(f"Variant count for {variant}: {variant_count}")
        snv_files = variant_df
      # Perform frequency count on insertion file
      elif variant == "insertions":
        logger.debug(f"Variant file of type {variant} selected")
        insertion_count = str((variant_df["Frequency"] >= variant_threshold).sum())
        logger.debug(f"Variant count for {variant}: {insertion_count}")
        insertion_files = variant_df
      # Perform frequency count on deletions file
      elif variant == "deletions":
        logger.debug(f"Variant file of type {variant} selected")
        deletion_count = str((variant_df["Frequency"] >= variant_threshold).sum())
        logger.debug(f"Variant count for {variant}: {deletion_count}")
        deletion_files = variant_df
    except IndexError as e:
      logger.error(f"Error parsing variant file {variant}: {e}")
      
  return deletion_count, insertion_count, variant_count, snv_files, insertion_files, deletion_files

def create_mira_qc(
                  input_dir: Path,
                  segments: Dict[str, str], 
                  samplename: str,
                  irma_type: str,
                  subtype_dict: Dict[str, str],
                  logger: logging.Logger 
):

  # Create header for QC table and declare dataframes for variant information
  qc_df = pd.DataFrame(columns=[
    'Sample', 'Total Reads', 'Pass QC', 'Reads Mapped', 'Reference',
    '% Reference Covered', 'Median Coverage', 'Mean Coverage',
    'Count of Minor SNVs (AF >= 0.05)',
    'Count of Minor Insertions (AF >= 0.05)',
    'Count of Minor Deletions (AF >= 0.05)'
  ])

  all_snv_files = pd.DataFrame()
  all_deletion_files = pd.DataFrame()
  all_insertion_files = pd.DataFrame()

  for segment_idx in segments:
    segment = segments[segment_idx]
    # Assign subtype suffix for path building with subtype dictionary from subtype_selection
    subtype_suffix = f"_{subtype_dict[segment]}" if segment in ["HA", "NA"] and subtype_dict[segment] else ""
    logger.debug(f"Subtype suffix used for pathbuilding: {subtype_suffix}")
    # Obtain read metrics from READ_COUNTS.tsv
    total_reads, pass_qc_reads, reads_mapped = read_counts(input_dir, samplename, logger, segment=segment)
    # Obtain the reference sequence length
    reference_path = input_dir / f"intermediate/0-ITERATIVE-REFERENCES/R0-{irma_type}_{segment}{subtype_suffix}.ref"
    segment_ref_len = reference_length(samplename, logger, segment, reference_path)
    # Obtain coverage information 
    coverage_path = input_dir / f"tables/{irma_type}_{segment}{subtype_suffix}-coverage.txt"
    mean_cov, median_cov, segment_pct_ref_cov = coverage_parsing(samplename, logger, segment, segment_ref_len, coverage_path)
    
    # Obtain variant counts and create variant dataframes for concatenation
    variant_path = input_dir / f"tables/{irma_type}_{segment}{subtype_suffix}"
    logger.info(f"Variant path: {variant_path}")
    deletion_count, insertion_count, variant_count, segment_snv_files, segment_insertion_files, segment_deletion_files  = variant_parsing(logger, variant_path)

    # Concatenate segment variant files to overall variant files
    all_snv_files = pd.concat([all_snv_files, segment_snv_files], ignore_index=True)
    all_deletion_files = pd.concat([all_deletion_files, segment_deletion_files], ignore_index=True)
    all_insertion_files = pd.concat([all_insertion_files, segment_insertion_files], ignore_index=True)

    # Append the row to the DataFrame
    row = {
      'Sample': samplename,
      'Total Reads': str(total_reads),
      'Pass QC': str(pass_qc_reads),
      'Reads Mapped': str(reads_mapped),
      'Reference': segment + subtype_suffix,
      '% Reference Covered': segment_pct_ref_cov,
      'Median Coverage': str(median_cov),
      'Mean Coverage': str(mean_cov),
      'Count of Minor SNVs (AF >= 0.05)': variant_count,    
      'Count of Minor Insertions (AF >= 0.05)': insertion_count,
      'Count of Minor Deletions (AF >= 0.05)': deletion_count
    }
    # Concatenate row to the QC dataframe
    qc_df = pd.concat([qc_df, pd.DataFrame([row])], ignore_index=True)
    logger.debug(f"Row for segment {segment} added.")

  # Write concatenated variant files
  try:
    all_snv_files.to_csv(f"{samplename}/tables/{samplename}_irma_all_variants.tsv", sep="\t", index=False)
    all_deletion_files.to_csv(f"{samplename}/tables/{samplename}_irma_all_deletions.tsv", sep="\t", index=False)
    all_insertion_files.to_csv(f"{samplename}/tables/{samplename}_irma_all_insertions.tsv", sep="\t", index=False)
    logger.debug("Concatenated variant files written.")
  except IOError as e:
    logger.error(f"Error writing variant files: {e}")

  # Fill empty entries with N/A and write QC summary file
  qc_df.replace(['', None], "N/A", inplace=True)
  qc_df.to_csv(f"{samplename}/{samplename}_irma_qc_summary.tsv", sep="\t", index=False)
  logger.debug("QC summary written.")

def main():
  parser = argparse.ArgumentParser(
      description="A helper script to parse IRMA output files and create various consensus FASTA files",
      formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )

  parser.add_argument("-d", "--input_dir", type=Path, required=True, help="Output directory of IRMA used as input for helper script")
  parser.add_argument("-t", "--irma_type", type=str, required=True, help="IRMA Flu Type; A or B")
  parser.add_argument("-s", "--samplename", type=str, required=True, help="Samplename used to build filepaths")

  # Logging configurations
  parser.add_argument("--log", default="irma_helper.log", help="Log file")
  parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
  args = parser.parse_args()
  
  logger = setup_logger(args.log, args.verbose)

  logger.info(f"Input Directory: {args.input_dir}")
  logger.info(f"Input IRMA Type: {args.irma_type}")
  logger.info(f"Samplename: {args.samplename}")

  logger.info("Running subtype selection")
  segments_dict, subtype = subtype_selection(args.input_dir, args.samplename, args.irma_type, logger)
  logger.info("Running consensus creation")
  consensus_creation(args.input_dir, args.samplename, segments_dict, logger)
  logger.info("Creating MIRA QC output")
  create_mira_qc(args.input_dir, segments_dict, args.samplename, args.irma_type, subtype, logger)

if __name__ == "__main__":
  main()