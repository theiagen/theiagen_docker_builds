#!/usr/bin/env python3
import sys
import argparse
import logging
import os
import glob
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from typing import Dict
"""
A helper script to facilitate QC summary creation and FASTA concatenation.
"""

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

def segment_selection(
  samplename: str,
  irma_type: str,
  logger: logging.Logger
):
  # flu segments from largest to smallest
  # Thank you Molly H. for this code block!
  # declare associative arrays for segment numbers
  full_type = ""
  subtype_notes = ""

  if irma_type == "Type_A":
    logger.info(f"Type reported as: {irma_type}")
    flu_segments = { "1": "PB2", "2": "PB1", "3": "PA", "4": "HA", "5": "NP", "6": "NA", "7": "MP", "8": "NS"}
    # Check for any subtypes within Type A segments, relying on structure of A_{segment}_{subtype}.fasta
    subtype_list = sorted(glob.glob(f"{samplename}/*_*_*.fasta"))
    logger.debug(f"List of files that contain subtype: {subtype_list}")

    if subtype_list:
      for file in subtype_list:
        subtype = Path(file).stem.split("_")[-1]
        full_type += subtype
      logger.info(f"Subtype reported as: {full_type}")
    else:
      full_type = "No subtype predicted by IRMA"
      logger.info(f"No subtype predicted by IRMA")

  elif irma_type == "Type_B":
    logger.info(f"Type reported as: {irma_type}")
    flu_segments: Dict[str, str] = { "2": "PB1", "1": "PB2", "3": "PA", "4": "HA", "5": "NP", "6": "NA", "7": "MP", "8": "NS"}
    subtype_notes = "IRMA does not differentiate Victoria and Yamagata Flu B lineages. See abricate_flu_subtype output column"
    logger.info(f"IRMA does not differentiate Victoria and Yamagata Flu B lineages. See abricate_flu_subtype output column")
  
  # Write subtype and subtype notes to file for WDL interpretation
  try:
    logger.debug("Writing IRMA_SUBTYPE to file.")
    with open("IRMA_SUBTYPE.txt", "w") as irma_subtype:
      irma_subtype.writelines(full_type)

    logger.info("Writing IRMA_SUBTYPE_NOTES to file.")
    with open("IRMA_SUBTYPE_NOTES.txt", "w") as subtype_notes_file:
      subtype_notes_file.writelines(subtype_notes)
  except Exception as e:
    logger.error(f"Error writing IRMA_SUBTYPE.txt file: {e}")

  return flu_segments, subtype

def consensus_creation(
  samplename: str,
  flu_segments: Dict[str, str],
  logger: logging.Logger
):
  
  consensus_array = []
  padded_consensus_array = []
  concatenated_seq = ""

  for segment in flu_segments:
    logger.debug(f"Parsing {flu_segments[segment]} in output files")
    segment_file = glob.glob(f"{samplename}/amended_consensus/*_{segment}.fa")
    if len(segment_file) > 0: 
      record = list(SeqIO.parse(segment_file[0], "fasta"))[0]
      record.id = record.id.replace(segment, flu_segments[segment])
      record.description = ""

      concatenated_seq += str(record.seq)
      consensus_array.append(record)

      SeqIO.write(record, segment_file[0], "fasta-2line")
      os.rename(segment_file[0], f"{samplename}/amended_consensus/{samplename}_{flu_segments[segment]}.fasta")
      logger.debug(f"Amended consensus FASTA for {samplename}_{flu_segments[segment]} written and renamed")

      # Replace . for N and remove -
      record.seq = Seq(str(record.seq).replace('.', 'N').replace('-', ''))
      SeqIO.write(record, f"padded_assemblies/{samplename}_{flu_segments[segment]}.pad.fasta", "fasta-2line")
      logger.debug(f"Padded fasta written for: {segment_file[0]}")
      
      padded_consensus_array.append(record)
    else:
      logger.warning(f"No file found for segment: {flu_segments[segment]}")

  try:
    SeqIO.write(consensus_array, f"{samplename}/amended_consensus/{samplename}.irma.consensus.fasta", "fasta-2line")
    SeqIO.write(padded_consensus_array, f"padded_assemblies/{samplename}.irma.consensus.pad.fasta", "fasta-2line")
    SeqIO.write(SeqRecord(Seq(concatenated_seq), f"{samplename}_irma_concatenated", description=""), f"{samplename}/amended_consensus/{samplename}.irma.consensus.concatenated.fasta", "fasta-2line")
    SeqIO.write(SeqRecord(Seq(str(concatenated_seq).replace('.', 'N').replace('-', '')), f"{samplename}_irma_concatenated_padded", description=""), f"padded_assemblies/{samplename}.irma.consensus.concatenated.pad.fasta", "fasta-2line")
  except Exception as e:
    logger.error(f"Error writing FASTA files: {e}")

def create_mira_qc(segments: Dict[str, str], 
                  samplename: str,
                  subtype: str,
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
  var_types = ["variants", "insertions", "deletions"]
  variant_threshold = 0.05
  snv_files = pd.DataFrame()
  deletion_files = pd.DataFrame()
  insertion_files = pd.DataFrame()

  # Read in READ_COUNTS.tsv and pull read metrics
  try:
    read_counts = pd.read_csv(f"{samplename}/tables/READ_COUNTS.tsv", sep='\t', index_col=0)
    logger.debug(f"READ_COUNTS.tsv file found for {samplename}.")
    total_reads = read_counts.loc["1-initial", "Reads"]
    pass_qc_reads = read_counts.loc["2-passQC", "Reads"]
    logger.info(f"Total Reads found to be: {total_reads} \n Passing Reads found to be: {pass_qc_reads}")
  except:
    logger.warning("WARNING: READ_COUNTS.tsv file not found for {samplename}. Cannot extract read counts for QC summary.")

  for segment in segments:
    variant_count = 0
    insertion_count = 0
    deletion_count = 0
    # Obtain mapped reads per segment
    if not read_counts.empty:
      logger.debug("READ_COUNTS.tsv present, parsing mapped reads")
      reads_mapped = read_counts.loc[read_counts.index.str.contains(segments[segment]), 'Reads'].values[0]
      logger.info(f"Mapped reads found to be: {reads_mapped}")

    # Load reference seq to obtain length of seq
    try:
      logger.debug(glob.glob(f"{samplename}/intermediate/0-ITERATIVE-REFERENCES/R0-*_{segments[segment]}*.ref")[0])
      ref_record = list(SeqIO.parse(glob.glob(f"{samplename}/intermediate/0-ITERATIVE-REFERENCES/R0-*_{segments[segment]}*.ref")[0], "fasta"))[0]
      segment_ref_len=str(len(ref_record))
      logger.debug(f"Segment reference found with length {segment_ref_len} for {samplename}:{segments[segment]}")
    except:
      logger.warning(f"WARNING: No reference file found for segment {segments[segment]} for {samplename}")
      segment_ref_len = "N/A"

    # Load coverage file for parsing
    try:
      coverages = pd.read_csv(glob.glob(f"{samplename}/tables/*_{segments[segment]}*-coverage.txt")[0], sep='\t')
      logger.debug(f"Coverage file found for {samplename}:{segments[segment]}")

      # Check for segment length presense 
      if segment_ref_len != "N/A" and int(segment_ref_len) > 0:
        mapped_bases = 0
        # Mapped bases
        for base in coverages['Consensus']:
          if base not in ['-','N', 'a', 't', 'c', 'g']:
            mapped_bases += 1
        # Calculate percentage rounded to second decimal place
        if mapped_bases > 0:
          segment_pct_ref_cov = round((mapped_bases / int(segment_ref_len) * 100), 2)
        # Obtain mean and median coverage from coverage file
        mean_cov = round(coverages['Coverage Depth'].mean(), 2)
        median_cov = round(coverages['Coverage Depth'].median(), 2)
        logger.info(f"Mean Coverage and Median Coverage calculated as {mean_cov} and {median_cov} respectively.")
    except Exception as e:
      logger.warning(f"WARNING: segment_coverage_file is not set. Cannot calculate coverage statistics for segment {segments[segment]} for {samplename}.")
      mean_cov = "N/A"
      median_cov = "N/A"
      segment_pct_ref_cov = None

    # Load variant tables
    try: 
      # Using var_types reference dictionary obtain variant frequency for each type per segment
      # and append to concatenated variant files
      for type in var_types:
        logger.debug(f"Checking type {type}")

        try:
          variant_file = glob.glob(f"{samplename}/tables/*_{segments[segment]}*-{type}.txt")[0]
        except:
          logger.warning(f"No files found for type {type}")
          continue

        logger.debug(f"Variant file selected: {variant_file}")
        variant_df = pd.read_csv(variant_file, sep='\t')
        logger.debug("Variant file loaded as dataframe")
        # Perform frequency count on variants file
        if type == "variants":
          logger.debug(f"Variant file of type {type} selected")
          variant_count = (variant_df["Minority_Frequency"] > variant_threshold).sum()
          logger.debug(f"Variant count for {type}: {variant_count}")
          snv_files = pd.concat([snv_files, variant_df], ignore_index=True)
        # Perform frequency count on insertion file
        elif type == "insertions":
          logger.debug(f"Variant file of type {type} selected")
          insertion_count = (variant_df["Frequency"] > variant_threshold).sum()
          logger.debug(f"Variant count for {type}: {insertion_count}")
          insertion_files = pd.concat([insertion_files, variant_df], ignore_index=True)
        # Perform frequency count on deletions file
        elif type == "deletions":
          logger.debug(f"Variant file of type {type} selected")
          deletion_count = (variant_df["Frequency"] > variant_threshold).sum()
          logger.debug(f"Variant count for {type}: {deletion_count}")
          deletion_files = pd.concat([deletion_files, variant_df], ignore_index=True)
    except Exception as e:
      logger.error(e)
          
    # Append the row to the DataFrame
    row = {
      'Sample': samplename,
      'Total Reads': str(total_reads),
      'Pass QC': str(pass_qc_reads),
      'Reads Mapped': str(reads_mapped),
      'Reference': segments[segment] + "_" + subtype,
      '% Reference Covered': str(segment_pct_ref_cov),
      'Median Coverage': str(median_cov),
      'Mean Coverage': str(mean_cov),
      'Count of Minor SNVs (AF >= 0.05)': str(variant_count),    
      'Count of Minor Insertions (AF >= 0.05)': str(insertion_count),
      'Count of Minor Deletions (AF >= 0.05)': str(deletion_count)
    }
    # Concatenate row to the QC dataframe
    qc_df = pd.concat([qc_df, pd.DataFrame([row])], ignore_index=True)
    logger.debug(f"Row for segment {segments[segment]} added.")

  # Write concatenated variant files
  try:
    snv_files.to_csv(f"{samplename}/tables/{samplename}_irma_all_variants.tsv", sep="\t", index=False)
    deletion_files.to_csv(f"{samplename}/tables/{samplename}_irma_all_deletions.tsv", sep="\t", index=False)
    insertion_files.to_csv(f"{samplename}/tables/{samplename}_irma_all_insertions.tsv", sep="\t", index=False)
    logger.debug("Concatenated variant files written.")
  except Exception as e:
    logger.warning(e)

  # Fill empty entries with N/A and write QC summary file
  qc_df.fillna("N/A", inplace=True)
  qc_df.to_csv(f"{samplename}/{samplename}_irma_qc_summary.tsv", sep="\t", index=False)
  logger.debug("QC summary written.")

def main():
  parser = argparse.ArgumentParser(
      description="A helper script to increase readability and modularity of IRMA WDL task",
      formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )

  parser.add_argument("-t", "--type", type=str, required=True, help="IRMA Flu Type; Type_A or Type_B")
  parser.add_argument("-s", "--samplename", type=str, required=True, help="Samplename used to identify output directory and filenames")
  
  # Logging configurations
  parser.add_argument("--log", default="irma_helper.log", help="Log file")
  parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
  args = parser.parse_args()
  
  logger = setup_logger(args.log, args.verbose)

  segments_dict, subtype = segment_selection(args.samplename, args.type, logger)
  consensus_creation(args.samplename, segments_dict, logger)
  create_mira_qc(segments_dict, args.samplename, subtype, logger)

if __name__ == "__main__":
  main()