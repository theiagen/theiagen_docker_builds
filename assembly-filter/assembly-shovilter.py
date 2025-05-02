#!/usr/bin/env python3
import sys
import re
import argparse
import logging
import gzip
from typing import List, Tuple, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

"""
Core contig filtering script based on Shovill's filtering logic.
Implements length filtering, coverage filtering, and homopolymer filtering.
Allows for filtering options to be skipped.
"""

def setup_logger(log_file: str, verbose: bool) -> logging.Logger:
    """Set up to track in the contig filtering task."""
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


def extract_coverage_from_header(header: str) -> Optional[float]:
    """
    Extract coverage information from contig header.
    
    Handles formats from various assemblers per Shovill's logic:
    - SPAdes: >NODE_1_length_114969_cov_29.8803_pilon
    - MEGAHIT: >k127_7 flag=0 multi=23.8788 len=22
    - SKESA: >Contig_2_53.4039
    
    Args:
        header: The contig header string
        
    Returns:
        Float coverage value or None if not found
    """
    match = re.search(r'(multi=|cov_|Contig_\d+_)(\d+(\.\d+)?)', header)
    return float(match.group(2)) if match else None


def filter_contigs_by_length(
    records: List[SeqRecord], 
    min_length: int,
    logger: logging.Logger
) -> Tuple[List[SeqRecord], int]:
    """
    Filter contigs by minimum length.
    
    Args:
        records: List of SeqRecord objects
        min_length: Minimum contig length
        logger: Logger object
        
    Returns:
        Tuple of (filtered records, number of filtered contigs)
    """
    filtered_records = []
    filtered_count = 0
    
    for record in records:
        seq_len = len(record.seq)
        
        if seq_len < min_length:
            logger.info(f"Filtered contig {record.id}: Too short ({seq_len} < {min_length})")
            filtered_count += 1
            continue
        
        filtered_records.append(record)
    
    return filtered_records, filtered_count


def filter_contigs_by_coverage(
    records: List[SeqRecord], 
    min_coverage: float,
    logger: logging.Logger
) -> Tuple[List[SeqRecord], int]:
    """
    Filter contigs by minimum coverage.
    
    Args:
        records: List of SeqRecord objects
        min_coverage: Minimum contig coverage
        logger: Logger object
        
    Returns:
        Tuple of (filtered records, number of filtered contigs)
    """
    filtered_records = []
    filtered_count = 0
    
    for record in records:
        coverage = extract_coverage_from_header(record.description)
        
        if coverage is not None and coverage < min_coverage:
            logger.info(f"Filtered contig {record.id}: Low coverage ({coverage:.1f} < {min_coverage})")
            filtered_count += 1
            continue
        elif coverage is None:
            logger.warning(f"Contig {record.id} has no coverage information, skipping coverage filter")
        
        filtered_records.append(record)
    
    return filtered_records, filtered_count


def filter_contig_homopolymers(
    records: List[SeqRecord],
    logger: logging.Logger
) -> Tuple[List[SeqRecord], int]:
    """
    Filter out homopolymer contigs from biopython SeqRecord objects.
    
    Args:
        records: List of SeqRecord objects
        logger: Logger object
        
    Returns:
        Tuple of (filtered records (no homopolymer), number of filtered contigs)
    """
    filtered_records = []
    filtered_count = 0
    
    for record in records:
        seq_str = str(record.seq).upper()
        # Check for homopolymers and filter them out
        if len(set(seq_str)) == 1:
            # If all bases are the same, it's a homopolymer, we can grab the first base as its the same
            # for all bases in the contig
            base = seq_str[0]
            logger.info(f"Filtered contig {record.id}: Homopolymer ({len(seq_str)} bases of {base})")
            filtered_count += 1
            # Skip this record and do not add it to the filtered records
            continue
        
        filtered_records.append(record)
    
    return filtered_records, filtered_count


def main():
    parser = argparse.ArgumentParser(
        description="Filter contigs based on length, coverage, and homopolymers",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-i", "--input", required=True,
                       help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True,
                       help="Output FASTA file")
    parser.add_argument("-m", "--metrics", default="filtering_metrics.txt",
                       help="Metrics output file")
    # Optional arguments for core filtering options in case they need to be skipped
    parser.add_argument("--skip-length-filter", action="store_true",
                       help="Disable length filtering")
    parser.add_argument("--skip-coverage-filter", action="store_true",
                       help="Disable coverage filtering")
    parser.add_argument("--skip-homopolymer-filter", action="store_true",
                       help="Disable homopolymer filtering")
    # Parameters for contig filtering
    parser.add_argument("--minlen", type=int, default=500,
                       help="Minimum contig length")
    parser.add_argument("--mincov", type=float, default=2.0,
                       help="Minimum contig coverage")
    # Logging configurations
    parser.add_argument("--log", default="contig_filter.log",
                       help="Log file")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Verbose output")
    
    args = parser.parse_args()
    
    logger = setup_logger(args.log, args.verbose)
    
    # Read in the fasta file and process it with biopython
    try:
        logger.info(f"Reading input file: {args.input}")
        if args.input.endswith('.gz'):
            with gzip.open(args.input, 'rt') as handle:
                records = list(SeqIO.parse(handle, "fasta"))
        else:
            records = list(SeqIO.parse(args.input, "fasta"))
        total_contigs = len(records)
        total_bases = sum(len(record.seq) for record in records)
        logger.info(f"Read {total_contigs} contigs, {total_bases} bases")
        
        # Track metrics for reporting like is done in flye assembly
        length_filtered = 0
        coverage_filtered = 0
        homopolymer_filtered = 0
        filtered_records = records
        
        # We need to be able to skip the filtering steps if the user wants to, 
        # i.e. coverage filtering is not always needed, but most scenarios will want length filtering and homopolymer filtering
        if not args.skip_length_filter:
            logger.info(f"Filtering contigs by length (min {args.minlen} bp)")
            filtered_records, length_filtered = filter_contigs_by_length(
                filtered_records, args.minlen, logger
            )
            logger.info(f"Removed {length_filtered} contigs with length < {args.minlen} bp")
        
        if not args.skip_coverage_filter:
            logger.info(f"Filtering contigs by coverage (min {args.mincov}x)")
            filtered_records, coverage_filtered = filter_contigs_by_coverage(
                filtered_records, args.mincov, logger
            )
            logger.info(f"Removed {coverage_filtered} contigs with coverage < {args.mincov}x")
        
        if not args.skip_homopolymer_filter:
            logger.info("Filtering homopolymer contigs")
            filtered_records, homopolymer_filtered = filter_contig_homopolymers(
                filtered_records, logger
            )
            logger.info(f"Removed {homopolymer_filtered} homopolymer contigs")
        
        # We want ot write out the filtered records, like is done currently in the assemby filter task
        retained_contigs = len(filtered_records)
        retained_bases = sum(len(record.seq) for record in filtered_records)
        logger.info(f"Writing {retained_contigs} contigs to {args.output}")
        SeqIO.write(filtered_records, args.output, "fasta")
        
        # Write metrics -- like is done in the assembly filter task
        with open(args.metrics, "w") as metrics_file:
            metrics_file.write("Contig Filtering Metrics\n")
            metrics_file.write("========================\n")
            metrics_file.write(f"Total contigs before filtering: {total_contigs}\n")
            metrics_file.write(f"Total sequence length before filtering: {total_bases} bases\n")
            metrics_file.write(f"Total contigs after filtering: {retained_contigs}\n")
            metrics_file.write(f"Total sequence length after filtering: {retained_bases} bases\n")
            metrics_file.write(f"Contigs removed (short length): {length_filtered}\n")
            metrics_file.write(f"Contigs removed (low coverage): {coverage_filtered}\n")
            metrics_file.write(f"Contigs removed (homopolymers): {homopolymer_filtered}\n")
            metrics_file.write("========================\n")
            metrics_file.write("Filtering completed successfully.\n")
        
        logger.info("Filtering completed successfully")
        
    except Exception as exc:
        logger.error(f"Error: {exc}, exiting program", exc_info=args.verbose)
        sys.exit(1)

if __name__ == "__main__":
    main()