#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import os
import sys
import json
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(filename)s - %(levelname)s - %(message)s')

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def parse_args():
    parser = argparse.ArgumentParser(description='Prepare ENA reads metadata spreadsheet from Terra table')
    
    # Input/output arguments
    parser.add_argument('--input', required=True, help='Input TSV file exported from Terra')
    parser.add_argument('--output', required=True, help='Output TSV file for ENA submission')
    parser.add_argument('--file-paths', required=True, help='Output JSON file mapping original paths to local paths')
    parser.add_argument('--excluded', default="excluded_samples.tsv", help='Output file for excluded samples')
    
    # Sample selection
    parser.add_argument('--sample-id-column', required=True, help='Column name containing sample IDs')
    parser.add_argument('--samples', required=False, help='Comma-separated list of sample IDs to include (default: all)')
    
    # Column mappings if user wants to map Terra columns to ENA columns
    parser.add_argument('--column-mappings', help='TSV file containing mappings from Terra columns to ENA columns')
    
    # Column name configuration (used if column-mappings not provided)
    parser.add_argument('--read1-column', default='read1', help='Column name for read1/first fastq files (default: read1)')
    parser.add_argument('--read2-column', default='read2', help='Column name for read2/second fastq files (default: read2)')
    parser.add_argument('--bam-column', default='bam_file', help='Column name for BAM files (default: bam_file)')
    parser.add_argument('--cram-column', default='cram_file', help='Column name for CRAM files (default: cram_file)')
    
    # Pass in if not in the table, study_accession always required -- similar to how terra2_ncbi works
    parser.add_argument('--study-accession', required=True, help='ENA study accession')
    parser.add_argument('--experiment-name', help='Default experiment name template')
    parser.add_argument('--library-strategy', default='WGS', help='Default library strategy (default: WGS)')
    parser.add_argument('--library-source', default='GENOMIC', help='Default library source (default: GENOMIC)')
    parser.add_argument('--library-selection', default='RANDOM', help='Default library selection (default: RANDOM)')
    parser.add_argument('--platform', default='ILLUMINA', help='Default sequencing platform (default: ILLUMINA)')
    parser.add_argument('--instrument', help='Default sequencing instrument')
    
    # Directory for storing data files
    parser.add_argument('--data-dir', default='data', help='Directory to store data files (default: data)')
    
    # Allow continuing even if some samples have missing metadata
    # Default is to fail if any samples are missing required metadata, unless this flag is set
    parser.add_argument('--allow-missing', action='store_true', 
                       help='Allow processing to continue even if some samples have missing required metadata')
    
    return parser.parse_args()

def check_required_metadata(table, required_metadata, sample_id_column):
    """
    Check if any rows are missing required metadata
    
    Args:
        table: The table to check
        required_metadata: List of required column names
        sample_id_column: Name of the sample ID column
        
    Returns:
        tuple: (missing_data_exists, excluded_samples_df)
    """
    # Replace blank cells with NaNs
    table.replace(r'^\s+$', np.nan, regex=True)
    
    # Find rows with missing required data
    excluded_samples = table[table[required_metadata].isna().any(axis=1)].copy()
    
    if len(excluded_samples) > 0:
        # Format the excluded samples for reporting
        excluded_samples.set_index(sample_id_column, inplace=True)
        excluded_samples = excluded_samples[excluded_samples.columns.intersection(required_metadata)]
        excluded_samples = excluded_samples.loc[:, excluded_samples.isna().any()]
        return True, excluded_samples
    
    return False, pd.DataFrame()

def main():
    args = parse_args()
    
    # Load the input table
    try:
        table = pd.read_csv(args.input, sep='\t', header=0, dtype=str)
        logger.info(f"Loaded table with {len(table)} rows and {len(table.columns)} columns")
    except Exception as e:
        logger.info(f"Error loading input table: {e}")
        sys.exit(1)
    
    # Load column mappings if provided
    column_mappings = {}
    if args.column_mappings:
        try:
            mappings_df = pd.read_csv(args.column_mappings, sep='\t', header=0)
            if 'terra_column' in mappings_df.columns and 'ena_column' in mappings_df.columns:
                for _, row in mappings_df.iterrows():
                    if pd.notna(row['terra_column']) and pd.notna(row['ena_column']):
                        column_mappings[row['terra_column']] = row['ena_column']
                logger.info(f"Loaded {len(column_mappings)} column mappings")
            else:
                logger.info("Warning: Column mappings file must have 'terra_column' and 'ena_column' columns")
        except Exception as e:
            logger.info(f"Error loading column mappings: {e}")
    
    # Set up the file column mappings based on input args or column mappings file
    read1_column = args.read1_column
    read2_column = args.read2_column
    bam_column = args.bam_column
    cram_column = args.cram_column
    
    # Override with mappings from file if provided
    if column_mappings:
        # Look for read file columns in mappings
        for terra_col, ena_col in column_mappings.items():
            if ena_col == 'uploaded file 1' and terra_col in table.columns:
                read1_column = terra_col
            elif ena_col == 'uploaded file 2' and terra_col in table.columns:
                read2_column = terra_col
            elif ena_col == 'bam_file' and terra_col in table.columns:
                bam_column = terra_col
            elif ena_col == 'cram_file' and terra_col in table.columns:
                cram_column = terra_col
    
    # Filter to specified samples if provided
    if args.samples:
        sample_list = args.samples.split(',')
        table = table[table[args.sample_id_column].isin(sample_list)]
        logger.info(f"Filtered to {len(table)} samples")
    
    # Required metadata fields for ENA submission - raw reads
    #https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html
    required_metadata = [
        "sample_accession",
        "experiment_name",
        "sequencing_platform",
        "sequencing_instrument",
        "library_source",
        "library_selection",
        "library_strategy"
    ]
    
    # Need at least one read file
    read_file_column = None
    if read1_column in table.columns:
        required_metadata.append(read1_column)
        read_file_column = read1_column
    elif bam_column in table.columns:
        required_metadata.append(bam_column)
        read_file_column = bam_column
    elif cram_column in table.columns:
        required_metadata.append(cram_column)
        read_file_column = cram_column
    
    if not read_file_column:
        logger.error(f"Error: No read file columns found in the table. At least one of {read1_column}, {bam_column}, or {cram_column} must be present.")
        sys.exit(1)
    else:
        logger.info(f"Using {read_file_column} as the primary read file column")
    
    # Check for missing required metadata and fail if any is found (unless --allow-missing is set)
    missing_data, excluded_samples = check_required_metadata(table, required_metadata, args.sample_id_column)
    
    if missing_data:
        # Write out the excluded samples report
        with open(args.excluded, "w") as exclusions:
            exclusions.write(f"Samples excluded for missing required metadata for reads submission (will have empty values in indicated columns):\n")
        
        excluded_samples.to_csv(args.excluded, mode='a', sep='\t')
        logger.error(f"{len(excluded_samples)} samples are missing required metadata:")
        logger.error(excluded_samples)
        
        if not args.allow_missing:
            logger.error("\nSubmission aborted. All samples must have required metadata.")
            logger.error("Required metadata fields are: " + ", ".join(required_metadata))
            logger.error("\nYou can override this check with --allow-missing but only samples with complete metadata will be submitted.")
            sys.exit(1)
        else:
            logger.info("\nContinuing with remaining samples as --allow-missing was specified.")
            # Remove rows with missing metadata
            table = table.dropna(subset=required_metadata, axis=0, how='any')
            logger.info(f"{len(table)} samples remain after removing those with missing required metadata")
    else:
        # If no missing data, write an empty excluded samples report
        with open(args.excluded, "w") as exclusions:
            exclusions.write("No samples were excluded for missing required metadata\n")
        logger.info("All samples have the required metadata")
    
    # Track file paths and their locations
    file_paths = {}
    
    # Create a new dataframe with exactly the columns expected by the ENA bulk submission tool
    columns = [
        "study_accession", "sample_accession", "experiment_name", "sequencing_platform",
        "sequencing_instrument", "library_name", "library_source", "library_selection",
        "library_strategy", "library_description", "insert_size", "uploaded file 1", "uploaded file 2"
    ]
    
    ena_spreadsheet = pd.DataFrame(columns=columns)
    
    # Create a mapping of ENA columns to Terra columns based on column_mappings
    ena_to_terra_map = {}
    if column_mappings:
        for terra_col, ena_col in column_mappings.items():
            ena_to_terra_map[ena_col] = terra_col
    
    # Process each row in the table to create the ENA spreadsheet
    for index, row in table.iterrows():
        new_row = {}
        
        # Study accession - use input parameter
        new_row["study_accession"] = args.study_accession
        
        # Sample accession - required
        sample_accession_col = ena_to_terra_map.get("sample_accession", "sample_accession")
        new_row["sample_accession"] = row[sample_accession_col]
        
        # Process other columns using the mapping
        for ena_col in columns:
            if ena_col in ["study_accession", "sample_accession", "uploaded file 1", "uploaded file 2"]:
                # These are handled separately
                continue
            
            # First check to see if there is a mapping for this column in the column mappings file
            # Verify that the mapped column exists in the row and has a value
            # If not, use defaults or generate values
            if ena_col in ena_to_terra_map and ena_to_terra_map[ena_col] in row and pd.notna(row[ena_to_terra_map[ena_col]]):
                # Get value from the mapped Terra column
                new_row[ena_col] = row[ena_to_terra_map[ena_col]]
            else:
                # Use defaults or generate values
                if ena_col == "experiment_name":
                    if "experiment_name" in row and pd.notna(row["experiment_name"]):
                        new_row[ena_col] = row["experiment_name"]
                    else:
                        new_row[ena_col] = args.experiment_name
                
                # Don't need sequencing_platform if we have instrument_model, per docs
                elif ena_col == "sequencing_platform":
                    if "sequencing_platform" in row and pd.notna(row["sequencing_platform"]):
                        new_row[ena_col] = row["sequencing_platform"]
                    elif args.platform:
                        new_row[ena_col] = args.platform
                    elif "sequencing_instrument" in row and pd.notna(row["instrument_model"]):
                        new_row[ena_col] = ""
                
                elif ena_col == "sequencing_instrument":
                    if "sequencing_instrument" in row and pd.notna(row["sequencing_instrument"]):
                        new_row[ena_col] = row["sequencing_instrument"]
                    else:
                        new_row[ena_col] = args.instrument
                
                # Library name is optional
                elif ena_col == "library_name":
                    if "library_name" in row and pd.notna(row["library_name"]):
                        new_row[ena_col] = row["library_name"]
                    else:
                        new_row[ena_col] = ""
                
                elif ena_col == "library_source":
                    if "library_source" in row and pd.notna(row["library_source"]):
                        new_row[ena_col] = row["library_source"]
                    else:
                        new_row[ena_col] = args.library_source
                
                elif ena_col == "library_selection":
                    if "library_selection" in row and pd.notna(row["library_selection"]):
                        new_row[ena_col] = row["library_selection"]
                    else:
                        new_row[ena_col] = args.library_selection
                
                elif ena_col == "library_strategy":
                    if "library_strategy" in row and pd.notna(row["library_strategy"]):
                        new_row[ena_col] = row["library_strategy"]
                    else:
                        new_row[ena_col] = args.library_strategy
                
                # This is optional
                elif ena_col == "library_description":
                    if "description" in row and pd.notna(row["description"]):
                        new_row[ena_col] = row["description"]
                    else:
                        new_row[ena_col] = ""
                
                # This is optional
                elif ena_col == "insert_size":
                    if "insert_size" in row and pd.notna(row["insert_size"]):
                        new_row[ena_col] = row["insert_size"]
                    else:
                        new_row[ena_col] = ""
        
        # Read files with path mapping
        if read1_column in row and pd.notna(row[read1_column]):
            original_path1 = row[read1_column]
            
            # Create a local filename from the original path
            filename1 = os.path.basename(original_path1)
            local_path1 = f"{args.data_dir}/{filename1}"
            
            # Store mapping of original to local path
            file_paths[original_path1] = local_path1
            
            # Use local path in spreadsheet
            new_row["uploaded file 1"] = local_path1
            
            if read2_column in row and pd.notna(row[read2_column]):
                original_path2 = row[read2_column]
                
                # Create a local filename from the original path
                filename2 = os.path.basename(original_path2)
                local_path2 = f"{args.data_dir}/{filename2}"
                
                # Store mapping of original to local path
                file_paths[original_path2] = local_path2
                
                # Use local path in spreadsheet
                new_row["uploaded file 2"] = local_path2
            else:
                new_row["uploaded file 2"] = ""
        elif bam_column in row and pd.notna(row[bam_column]):
            original_path = row[bam_column]
            
            # Create a local filename from the original path
            filename = os.path.basename(original_path)
            local_path = f"{args.data_dir}/{filename}"
            
            # Store mapping of original to local path
            file_paths[original_path] = local_path
            
            # Use local path in spreadsheet
            new_row["uploaded file 1"] = local_path
            new_row["uploaded file 2"] = ""
        elif cram_column in row and pd.notna(row[cram_column]):
            original_path = row[cram_column]
            
            # Create a local filename from the original path
            filename = os.path.basename(original_path)
            local_path = f"{args.data_dir}/{filename}"
            
            # Store mapping of original to local path
            file_paths[original_path] = local_path
            
            # Use local path in spreadsheet
            new_row["uploaded file 1"] = local_path
            new_row["uploaded file 2"] = ""
        
        # Add the row to the dataframe
        ena_spreadsheet = pd.concat([ena_spreadsheet, pd.DataFrame([new_row])], ignore_index=True)
    
    # Check if we have any samples left to submit
    if len(ena_spreadsheet) == 0:
        logger.info("No samples remain after filtering out those with missing metadata")
        sys.exit(1)
    
    # Write the ENA spreadsheet
    ena_spreadsheet.to_csv(args.output, sep="\t", index=False)
    logger.info(f"Created ENA reads spreadsheet with {len(ena_spreadsheet)} samples")
    
    # Write file paths mapping as JSON
    with open(args.file_paths, "w") as f:
        json.dump(file_paths, f, indent=2)
    logger.info(f"Listed {len(file_paths)} data files that need to be localized")

if __name__ == "__main__":
    main()