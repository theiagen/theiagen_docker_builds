#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import sys
import logging

# Small script, so using gsutil and cp directly instead of using a library
# This script is used to localize files from a remote path (gs://) to a local path
# The file paths are provided in a JSON file that maps the original paths to the local paths
# This script is used in the pipeline after the files have been prepared for submission to ENA

# Set up logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(filename)s - %(levelname)s - %(message)s')

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def parse_args():
    parser = argparse.ArgumentParser(description='Localize files for ENA submission')
    
    # File path from prepare_ena_data.py that contains the mapping of original paths to local paths which bulk_webincli will use
    parser.add_argument('--file-paths', required=True, help='JSON file mapping original paths to local paths')
    
    return parser.parse_args()

def copy_file(source_path, target_path):
    """Copy a file from source to target path"""
    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    
    # Check if source is a gs:// path
    if source_path.startswith('gs://'):
        logger.info(f"Downloading {source_path} to {target_path}")
        process = subprocess.run(['gsutil', 'cp', source_path, target_path], 
                                 capture_output=True, text=True)
        if process.returncode != 0:
            logger.error(f"Error downloading {source_path}: {process.stderr}")
            return False
    else:
        # Local file
        logger.info(f"Copying {source_path} to {target_path}")
        try:
            target_dir = os.path.dirname(target_path)
            if target_dir and not os.path.exists(target_dir):
                os.makedirs(target_dir)
                
            subprocess.run(['cp', source_path, target_path], check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Error copying {source_path}: {e}")
            return False
        except Exception as e:
            logger.error(f"Unexpected error copying {source_path}: {e}")
            return False
    
    return True

def main():
    args = parse_args()
    try:
        with open(args.file_paths, 'r') as f:
            file_paths = json.load(f)
        logging.info(f"Loaded {len(file_paths)} file paths to localize")
    except Exception as e:
        logger.info(f"Error loading file paths: {e}")
        sys.exit(1)
    
    # Track success and failures
    success_count = 0
    failure_count = 0
    
    # Copy each file
    for source_path, target_path in file_paths.items():
        if copy_file(source_path, target_path):
            success_count += 1
        else:
            failure_count += 1
    
    logger.info(f"Localization complete: {success_count} successful, {failure_count} failed")
    
    if failure_count > 0:
        logger.error(f"{failure_count} files failed to localize")
        sys.exit(1)

if __name__ == "__main__":
    main()