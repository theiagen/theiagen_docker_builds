#!/usr/bin/env python3

import argparse
import logging
import os
import pandas as pd
import requests
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from xml.dom import minidom
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(filename)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ENA API endpoints
ENDPOINT_PROD = {
    'dropbox': 'https://www.ebi.ac.uk/ena/submit/drop-box/submit/',
    'query': 'https://www.ebi.ac.uk/ena/portal/api/search?result=sample&query=(sample_alias={alias})'
}

ENDPOINT_TEST = {
    'dropbox': 'https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/',
    'query': 'https://wwwdev.ebi.ac.uk/ena/portal/api/search?result=sample&query=(sample_alias={alias})'
}

# NCBI API URL for taxonomy lookups
# If user provides a taxon_id we will use that, otherwise we will use the organism name to get the taxon_id
# Funny enough ENA depends on NCBI for taxon_id lookups...
NCBI_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{encoded_organism}"

# Checklist mappings for sample types
CHECKLIST_MAPPING = {
    "prokaryotic_pathogen": "ERC000028",
    "virus_pathogen": "ERC000033"
}

# Required fields for each checklist
CHECKLIST_REQUIREMENTS = {
    # ERC000028 - Prokaryotic pathogen minimal sample checklist
    "ERC000028": {
        "mandatory": [
            "collection date",
            "geographic location (country and/or sea)",
            "isolation_source", 
            "host scientific name",
            "host health state",
            "isolate"
        ],
        "recommended": [
            "lat_lon",
            "strain",
            "serovar"
        ]
    },
    # ERC000033 - Virus pathogen reporting standard checklist
    "ERC000033": {
        "mandatory": [
            "collection date",
            "geographic location (country and/or sea)",
            "host common name",
            "host scientific name",
            "host subject id",
            "host health state",
            "host sex",
            "isolate",
            "collector name", 
            "collecting institution"
        ],
        "recommended": [
            "geographic location (latitude)",
            "geographic location (longitude)",
            "geographic location (region and locality)",
            "sample capture status",
            "host behaviour",
            "host habitat",
            "host age",
            "isolation source host-associated",
            "isolation source non-host-associated",
            "virus identifier",
            "receipt date",
            "serotype (required for a seropositive sample)",
            "host disease outcome"
        ]
    }
}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Submit sample metadata to ENA using XML API')
    
    # Required arguments
    parser.add_argument('--metadata', required=True, help='Metadata spreadsheet (TSV format)')
    parser.add_argument('--output', required=True, help='Output directory for results')
    parser.add_argument('--username', required=True, help='ENA Webin username')
    parser.add_argument('--password', required=True, help='ENA Webin password')
    parser.add_argument('--study', required=True, help='ENA study accession')
    # In WDL this will come from entity_type from extract large tsv import terra-tool
    parser.add_argument('--sample-id-column', required=True, help='Column name containing sample IDs')
    parser.add_argument('--sample-type', required=True, 
                       choices=list(CHECKLIST_MAPPING.keys()), 
                       help='ENA sample type')
    
    # Optional arguments
    # We will want to clip off the read1 and read2 columns from the metadata file
    parser.add_argument('--read1-column', default='read1', help='Column name for read1/first fastq files (default: read1)')
    parser.add_argument('--read2-column', default='read2', help='Column name for read2/second fastq files (default: read2)')
    parser.add_argument('--center', help='Center name for submission')
    parser.add_argument('--test', action='store_true', help='Use test ENA server')
    parser.add_argument('--allow-missing', action='store_true', help='Allow missing metadata fields')
    parser.add_argument('--column-mappings', help='TSV file with column name mappings from Terra to ENA')
    parser.add_argument('--batch-size', type=int, default=100, 
                        help='Number of samples to submit in one batch (default: 100)')
    
    return parser.parse_args()

def load_column_mappings(mapping_file):
    """Load column name mappings from TSV file"""
    if not mapping_file or not os.path.exists(mapping_file):
        logger.warning(f"Column mapping file not found or not specified: {mapping_file}")
        return {}
    
    try:
        logger.info(f"Loading column mappings from {mapping_file}")
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        
        # Validate column names
        required_columns = ['terra_column', 'ena_column']
        if not all(col in mapping_df.columns for col in required_columns):
            logger.error(f"Column mapping file must have columns: {', '.join(required_columns)}")
            return {}
        
        # Create mapping dictionary
        mappings = {}
        for _, row in mapping_df.iterrows():
            terra_col = row['terra_column']
            ena_col = row['ena_column']
            
            if pd.notna(terra_col) and pd.notna(ena_col):
                mappings[terra_col] = ena_col
        
        logger.info(f"Loaded {len(mappings)} column mappings from {mapping_file}")
        return mappings
    
    except Exception as err:
        logger.error(f"Error loading column mappings file: {err}")
        return {}

def get_taxon_id_from_ncbi(organism):
    """Get NCBI Taxonomy ID by querying the NCBI Taxonomy API"""
    if not organism or str(organism).strip() in ('', 'NA', 'nan', 'None', 'NaN'):
        return None
    
    organism = str(organism).strip()
    
    try:
        # URL encode the organism name
        encoded_organism = requests.utils.quote(organism)
        url = NCBI_URL.format(encoded_organism=encoded_organism)
        
        # Make the request with timeout and retries
        for attempt in range(3):
            try:
                response = requests.get(url, timeout=10)
                break
            except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
                if attempt < 2:
                    logger.warning(f"NCBI API request timed out, retrying (attempt {attempt+1}/3)")
                time.sleep(1)  # Add a small delay between retries
        else:
            logger.warning(f"NCBI API request failed after 3 attempts for organism '{organism}'")
            return None
        
        # Check status code
        if response.status_code != 200:
            logger.warning(f"NCBI API returned status code {response.status_code} for organism '{organism}'")
            return None
        
        # Parse response
        data = response.json()
        
        # Get the taxon id from the first taxon id in the response
        if 'taxonomy_nodes' in data and len(data['taxonomy_nodes']) > 0:
            taxon_id = str(data['taxonomy_nodes'][0]['taxonomy']['tax_id'])
            logger.info(f"Found taxon ID {taxon_id} for organism '{organism}' via NCBI API")
            return taxon_id
        else:
            logger.warning(f"No taxonomy nodes found in NCBI API response for organism '{organism}'")
            return None
            
    except Exception as e:
        logger.warning(f"Error querying NCBI Taxonomy API for organism '{organism}': {e}")
        return None

def validate_metadata(metadata_df, checklist_id, sample_id_column, column_mappings=None, allow_missing=False):
    """Validate metadata against checklist requirements
    
    Returns:
        dict: Validation results with issues by sample
    """
    if checklist_id not in CHECKLIST_REQUIREMENTS:
        logger.warning(f"No requirements defined for checklist {checklist_id}")
        return {"valid": True, "issues": {}}
    
    requirements = CHECKLIST_REQUIREMENTS[checklist_id]
    validation_issues = {}
    
    # Default field mapping for common columns
    default_field_mapping = {
        'collection_date': 'collection date',
        'geo_loc_name': 'geographic location (country and/or sea)',
        'isolation_source': 'isolation_source',
        'host': 'host scientific name',
        'host_scientific_name': 'host scientific name',
        'host_common_name': 'host common name',
        'host_health_state': 'host health state',
        'host_sex': 'host sex',
        'host_subject_id': 'host subject id',
        'collector_name': 'collector name',
        'collecting_institution': 'collecting institution',
        'lat_lon': 'lat_lon',
        'latitude': 'geographic location (latitude)',
        'longitude': 'geographic location (longitude)',
        'region_locality': 'geographic location (region and locality)'
    }
    
    # Merge default with custom mappings
    field_mapping = default_field_mapping.copy()
    if column_mappings:
        field_mapping.update(column_mappings)
    
    for idx, row in metadata_df.iterrows():
        # I think these are all valid sample aliases, I guess if for somereason the sample_id_column is not in the row
        # then we can just use the idx to create a sample alias, but I think that is unlikely
        sample_alias = (row.get(sample_id_column) or
                       row.get('name') or 
                       row.get('submission_id') or 
                       row.get('sample_id') or 
                       f"sample_{idx}")
        
        # Create a combined row with both original and mapped column names
        combined_row = row.to_dict()
        for col, value in row.items():
            if col in field_mapping:
                combined_row[field_mapping[col]] = value
        
        # Check mandatory fields
        missing_mandatory = []
        for field in requirements.get("mandatory", []):
            # Check both the ENA field name and possible Terra column names
            found = False
            if field in combined_row and pd.notna(combined_row[field]) and str(combined_row[field]).strip() not in ('', 'NA', 'nan', 'None', 'NaN'):
                found = True
            
            if not found:
                # Try possible variations in column names
                for col, ena_field in field_mapping.items():
                    if ena_field == field and col in combined_row and pd.notna(combined_row[col]) and str(combined_row[col]).strip() not in ('', 'NA', 'nan', 'None', 'NaN'):
                        found = True
                        break
            
            if not found:
                missing_mandatory.append(field)
        
        # Check recommended fields
        missing_recommended = []
        for field in requirements.get("recommended", []):
            # Check both the ENA field name and possible Terra column names
            found = False
            if field in combined_row and pd.notna(combined_row[field]) and str(combined_row[field]).strip() not in ('', 'NA', 'nan', 'None', 'NaN'):
                found = True
            
            if not found:
                # Try possible variations in column names
                for col, ena_field in field_mapping.items():
                    if ena_field == field and col in combined_row and pd.notna(combined_row[col]) and str(combined_row[col]).strip() not in ('', 'NA', 'nan', 'None', 'NaN'):
                        found = True
                        break
            
            if not found:
                missing_recommended.append(field)
        
        if missing_mandatory:
            validation_issues[sample_alias] = {
                'missing_mandatory': missing_mandatory,
                'missing_recommended': missing_recommended
            }
            logger.warning(f"Sample {sample_alias} is missing mandatory fields: {', '.join(missing_mandatory)}")
        elif missing_recommended:
            logger.info(f"Sample {sample_alias} is missing recommended fields: {', '.join(missing_recommended)}")
    
    # Determine if validation passed
    is_valid = len(validation_issues) == 0 or allow_missing
    
    return {
        "valid": is_valid,
        "issues": validation_issues
    }

def prettify_xml(elem):
    """Return a pretty-printed XML string for the Element"""
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def add_sample_attribute(sample_attributes, tag, value, units=None):
    """Add a sample attribute to the XML"""
    if value is None or str(value).strip() in ('', 'NA', 'nan', 'None', 'NaN'):
        return False
    
    sample_attribute = ET.SubElement(sample_attributes, "SAMPLE_ATTRIBUTE")
    tag_elem = ET.SubElement(sample_attribute, "TAG")
    tag_elem.text = tag
    value_elem = ET.SubElement(sample_attribute, "VALUE")
    value_elem.text = str(value).strip()
    
    if units is not None and units != '':
        units_elem = ET.SubElement(sample_attribute, "UNITS")
        units_elem.text = units
    
    return True

def generate_sample_xml(metadata_df, checklist_id, sample_id_column, center_name=None, column_mappings=None):
    """Generate sample XML for ENA submission"""
    sample_set = ET.Element("SAMPLE_SET")
    
    # Field units for specific attributes
    field_units = {
        'depth': 'm',
        'elevation': 'm',
        'altitude': 'm',
        'temperature': 'celsius',
        'salinity': 'PSU',
        'geographic location (latitude)': 'DD',
        'geographic location (longitude)': 'DD',
        'host age': 'years'
    }
    
    # Default field mapping for common columns
    default_field_mapping = {
        'collection_date': 'collection date',
        'geo_loc_name': 'geographic location (country and/or sea)',
        'isolation_source': 'isolation_source',
        'host': 'host scientific name',
        'host_scientific_name': 'host scientific name',
        'host_common_name': 'host common name',
        'host_health_state': 'host health state',
        'host_sex': 'host sex',
        'host_disease': 'host disease',
        'host_disease_outcome': 'host disease outcome',
        'host_subject_id': 'host subject id',
        'host_age': 'host age',
        'host_behaviour': 'host behaviour',
        'host_habitat': 'host habitat',
        'isolation_source_host': 'isolation source host-associated',
        'environment_biome': 'environment (biome)',
        'environment_feature': 'environment (feature)',
        'environment_material': 'environment (material)',
        'collector_name': 'collector name',
        'collected_by': 'collector name',
        'collecting_institution': 'collecting institution',
        'receipt_date': 'receipt date',
        'sample_capture_status': 'sample capture status',
        'virus_identifier': 'virus identifier',
        'depth': 'depth',
        'elevation': 'elevation',
        'altitude': 'altitude',
        'temp': 'temperature',
        'salinity': 'salinity',
        'pH': 'pH',
        'habitat': 'habitat',
        'isolate': 'isolate',
        'strain': 'strain',
        'serovar': 'serovar',
        'serotype': 'serotype',
        'lab_host': 'lab_host',
        'latitude': 'geographic location (latitude)',
        'longitude': 'geographic location (longitude)',
        'region_locality': 'geographic location (region and locality)'
    }
    
    # Merge with custom mappings if provided
    field_mapping = default_field_mapping.copy()
    if column_mappings:
        # This will overwrite any default mappings with custom ones provided
        field_mapping.update(column_mappings)
    
    for idx, row in metadata_df.iterrows():
        # Get sample alias and title, same as in validate_metadata
        sample_alias = (row.get(sample_id_column) or
                        row.get('sample_accession') or 
                        row.get('name') or 
                        row.get('submission_id') or 
                        row.get('sample_id') or 
                        f"sample_{idx}")
        
        # Title should be passed in just like with NCBI
        title = (row.get('title'))
        if not title:
            logger.error(f"No title available for sample {sample_alias}")
            continue
        
        # Get organism and taxonomy information
        organism = row.get('organism')
        taxon_id = (row.get('taxon_id') or row.get('taxa_id') or row.get('taxonomy_id'))
        
        if not taxon_id and organism:
            taxon_id = get_taxon_id_from_ncbi(organism)
            if not taxon_id:
                logger.error(f"Failed to get taxon ID for '{organism}' from NCBI API")
                continue
        
        if not taxon_id:
            logger.error(f"No taxon ID available for sample {sample_alias}")
            continue
        
        # Create sample element
        sample = ET.SubElement(sample_set, "SAMPLE")
        sample.set("alias", str(sample_alias))
        
        if center_name:
            sample.set("center_name", center_name)
        
        # Add title
        title_elem = ET.SubElement(sample, "TITLE")
        title_elem.text = str(title)
        
        # Add taxonomy information
        sample_name = ET.SubElement(sample, "SAMPLE_NAME")
        taxon_id_elem = ET.SubElement(sample_name, "TAXON_ID")
        taxon_id_elem.text = str(taxon_id)
        
        if organism:
            scientific_name = ET.SubElement(sample_name, "SCIENTIFIC_NAME")
            scientific_name.text = str(organism)
        
        # Add description if available
        description = row.get('description') or row.get('library_description')
        if description and str(description).strip():
            description_elem = ET.SubElement(sample, "DESCRIPTION")
            description_elem.text = str(description)
        
        # Add sample attributes
        sample_attributes = ET.SubElement(sample, "SAMPLE_ATTRIBUTES")
        
        # First add ENA checklist
        add_sample_attribute(sample_attributes, "ENA-CHECKLIST", checklist_id)
        
        # Add ENA study link
        study_accession = row.get('study_accession')
        if study_accession:
            add_sample_attribute(sample_attributes, "ENA-STUDY", study_accession)
        
        # Once we get here we know all pontential fields are in the row and have been checked
        # so we can just loop through all the columns and add them to the sample attributes
        handled_fields = {'sample_accession', 'name', 'submission_id', 'title', 
                          'library_name', 'sample_title', 'organism', 'taxon_id',
                          'sample_type', 'attribute_package', 'study_accession', 
                          'description', 'library_description', 'sample_id'}
        
        # Add the dynamic sample_id_column to handled fields
        handled_fields.add(sample_id_column)
                          
        # Special handling for lat_lon
        lat_lon = row.get('lat_lon')
        if lat_lon and str(lat_lon).strip() not in ('', 'NA', 'nan', 'None', 'NaN'):
            try:
                # Try to parse lat,lon format assuming comma or space separator
                parts = str(lat_lon).strip().split(',')
                if len(parts) != 2:
                    parts = str(lat_lon).strip().split()
                
                if len(parts) == 2:
                    lat, lon = parts
                    # Add latitude
                    add_sample_attribute(sample_attributes, "geographic location (latitude)", lat.strip(), "DD")
                    # Add longitude
                    add_sample_attribute(sample_attributes, "geographic location (longitude)", lon.strip(), "DD")
                else:
                    # If parsing fails, just add as a single attribute
                    add_sample_attribute(sample_attributes, "geographic coordinates", str(lat_lon))
            except:
                # If parsing fails, just add as a single attribute
                add_sample_attribute(sample_attributes, "geographic coordinates", str(lat_lon))
            handled_fields.add('lat_lon')
        
        # Handle all mapped fields
        for col in row.index:
            # Skip already handled fields
            if col in handled_fields:
                continue
                
            value = row[col]
            if pd.isna(value) or str(value).strip() in ('', 'NA', 'nan', 'None', 'NaN'):
                continue
            
            # Use the mapped field name if available
            field_name = field_mapping.get(col, col)
            unit = field_units.get(field_name)
            
            add_sample_attribute(sample_attributes, field_name, value, unit)
    
    return sample_set

def generate_submission_xml(action="ADD"):
    """Generate submission XML for ENA
    Args:
        action (str): Action to perform (ADD, MODIFY, CANCEL)
    Returns:
        <?xml version="1.0" encoding="UTF-8"?>
        <SUBMISSION>
        <ACTIONS>
            <ACTION>
                <ADD/>
            </ACTION>
        </ACTIONS>
        </SUBMISSION>
    """
    submission = ET.Element("SUBMISSION")
    actions = ET.SubElement(submission, "ACTIONS")
    action_elem = ET.SubElement(actions, "ACTION")
    ET.SubElement(action_elem, action)
    return submission

def submit_to_ena(sample_xml, submission_xml, username, password, endpoint):
    """Submit XML to ENA and return response"""
    files = {
        "SAMPLE": ("sample.xml", ET.tostring(sample_xml, encoding="UTF-8", xml_declaration=True)),
        "SUBMISSION": ("submission.xml", ET.tostring(submission_xml, encoding="UTF-8", xml_declaration=True))
    }
    
    try:
        response = requests.post(
            endpoint['dropbox'],
            auth=(username, password),
            files=files
        )
        
        if response.status_code != 200:
            logger.error(f"HTTP error {response.status_code}: {response.text}")
            return None
        
        # Ensure we got XML response
        if "application/xml" not in response.headers.get("Content-Type", ""):
            logger.error(f"Expected XML response but got: {response.headers.get('Content-Type')}")
            logger.error(f"Response content: {response.text[:1000]}")
            return None
        
        return response.content
    except Exception as e:
        logger.error(f"Error submitting to ENA: {e}")
        return None

def process_receipt(receipt_xml):
    """Process receipt XML and extract accession numbers"""
    if receipt_xml is None:
        return None
    
    try:
        root = ET.fromstring(receipt_xml)
        success = root.get("success", "false").lower() == "true"
        
        results = []
        
        if success:
            # Extract sample accessions
            for sample in root.findall('.//SAMPLE'):
                sample_alias = sample.get('alias', '')
                ena_accession = sample.get('accession', '')
                status = sample.get('status', '')
                
                # Look for BioSample accession
                biosample_accession = ''
                for ext_id in sample.findall('.//EXT_ID'):
                    if ext_id.get('type') == 'biosample':
                        biosample_accession = ext_id.get('accession', '')
                
                results.append({
                    'sample_alias': sample_alias,
                    'ena_accession': ena_accession,
                    'biosample_accession': biosample_accession,
                    'status': status
                })
            
            return {
                'success': True,
                'results': results,
                'errors': []
            }
        else:
            # Extract error messages
            errors = []
            for message in root.findall('.//ERROR'):
                errors.append(message.text)
            
            return {
                'success': False,
                'results': [],
                'errors': errors
            }
    except Exception as e:
        logger.error(f"Error processing receipt XML: {e}")
        return {
            'success': False,
            'results': [],
            'errors': [str(e)]
        }

def batch_samples(metadata_df, batch_size):
    """Split samples into batches for submission"""
    num_samples = len(metadata_df)
    for i in range(0, num_samples, batch_size):
        yield metadata_df.iloc[i:min(i+batch_size, num_samples)]

def main():
    """Main function to process arguments and submit to ENA"""
    args = parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Set up logging to file
    file_handler = logging.FileHandler(os.path.join(args.output, 'submission.log'))
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    logger.info(f"Starting ENA submission, output directory: {args.output}")
    
    # Select endpoint based on test flag
    endpoint = ENDPOINT_TEST if args.test else ENDPOINT_PROD
    logger.info(f"Using {'TEST' if args.test else 'PRODUCTION'} ENA endpoint")
    
    # Load metadata spreadsheet
    try:
        metadata_df = pd.read_csv(args.metadata, sep='\t', dtype=str)
        metadata_df = metadata_df.replace(['NA', 'na', ''], pd.NA).fillna('')
        logger.info(f"Loaded metadata for {len(metadata_df)} samples")
        if args.sample_id_column not in metadata_df.columns:
            logger.error(f"Sample ID column '{args.sample_id_column}' not found in metadata file")
            logger.error(f"Available columns: {', '.join(metadata_df.columns)}")
            sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading metadata file: {e}")
        sys.exit(1)
    
    # Clip off readq1 and read2 columns, don't want to include in xml
    original_metadata_df = metadata_df.copy()
    metadata_drop_reads_df = original_metadata_df.drop(columns=[args.read1_column, args.read2_column])
    
    # Load column mappings if provided
    column_mappings = None
    if args.column_mappings:
        column_mappings = load_column_mappings(args.column_mappings)
    
    # Add study accession to metadata if not already there, just so we 
    # Grab everything we need from one place -- either way will be user provided on way or the other
    if 'study_accession' not in metadata_drop_reads_df.columns:
        metadata_df['study_accession'] = args.study
    
    # Get checklist ID for sample type
    checklist_id = CHECKLIST_MAPPING.get(args.sample_type)
    if not checklist_id:
        logger.error(f"Unknown sample type: {args.sample_type}")
        sys.exit(1)
    
    # Validate metadata against checklist
    validation_result = validate_metadata(
        metadata_drop_reads_df, 
        checklist_id,
        args.sample_id_column,
        column_mappings, 
        args.allow_missing
    )
    
    # Save validation results
    if validation_result['issues']:
        validation_file = os.path.join(args.output, 'validation_issues.json')
        with open(validation_file, 'w') as f:
            json.dump(validation_result['issues'], f, indent=2)
        
        logger.warning(f"Found metadata validation issues in {len(validation_result['issues'])} samples")
        logger.warning(f"See {validation_file} for details")
        
        if not validation_result['valid']:
            logger.error("Validation failed. Use --allow-missing to proceed with submission anyway.")
            with open(os.path.join(args.output, 'success.txt'), 'w') as f:
                f.write("false")
            return 1
    else:
        logger.info("All samples passed metadata validation!")
    
    # Process samples in batches
    all_results = []
    all_errors = []
    
    for batch_idx, batch_df in enumerate(batch_samples(metadata_drop_reads_df, args.batch_size)):
        logger.info(f"Processing batch {batch_idx+1}, samples {len(batch_df)}")
        
        # Generate sample XML
        sample_xml = generate_sample_xml(
            batch_df, 
            checklist_id,
            args.sample_id_column,
            args.center, 
            column_mappings
        )
        
        submission_xml = generate_submission_xml()
        
        # Save XML files
        xml_dir = os.path.join(args.output, 'xml')
        os.makedirs(xml_dir, exist_ok=True)
        
        batch_prefix = f"batch_{batch_idx+1}"
        sample_xml_path = os.path.join(xml_dir, f"{batch_prefix}_sample.xml")
        submission_xml_path = os.path.join(xml_dir, f"{batch_prefix}_submission.xml")
        
        with open(sample_xml_path, 'w') as f:
            f.write(prettify_xml(sample_xml))
        
        with open(submission_xml_path, 'w') as f:
            f.write(prettify_xml(submission_xml))
        
        # Submit to ENA
        logger.info(f"Submitting batch {batch_idx+1} to ENA")
        receipt_xml = submit_to_ena(sample_xml, submission_xml, args.username, args.password, endpoint)
        
        # Save receipt XML
        if receipt_xml is not None:
            receipt_xml_path = os.path.join(xml_dir, f"{batch_prefix}_receipt.xml")
            with open(receipt_xml_path, 'wb') as f:
                f.write(receipt_xml)
        
        # Process receipt
        receipt_result = process_receipt(receipt_xml)
        
        if receipt_result is None:
            logger.error(f"Failed to process receipt for batch {batch_idx+1}")
            all_errors.append(f"Failed to process receipt for batch {batch_idx+1}")
            continue
        
        # Save results
        if receipt_result['success']:
            all_results.extend(receipt_result['results'])
            logger.info(f"Successfully submitted {len(receipt_result['results'])} samples in batch {batch_idx+1}")
        else:
            all_errors.extend(receipt_result['errors'])
            logger.error(f"Errors in batch {batch_idx+1}: {', '.join(receipt_result['errors'])}")
        
        # Pause briefly between batches to avoid rate limiting just in case
        if batch_idx < len(list(batch_samples(metadata_df, args.batch_size))) - 1:
            time.sleep(2)
    
    # Write final results
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df.to_csv(os.path.join(args.output, 'accessions.tsv'), sep='\t', index=False)
        logger.info(f"Wrote {len(all_results)} accessions to accessions.tsv")
        
        # Create a mapping from sample_alias to accessions
        sample_accessions = {}
        for result in all_results:
            sample_alias = result['sample_alias']
            biosample_acc = result.get('biosample_accession', '')
            ena_acc = result.get('ena_accession', '')
            
            sample_accessions[sample_alias] = {
                'biosample_accession': biosample_acc,
                'ena_accession': ena_acc,
            }
        
        # Add accessions back to original metadata
        # Restore the original metadata with read columns
        updated_metadata = original_metadata_df.copy()
        
        # Update the accessions based on sample alias
        for idx, row in updated_metadata.iterrows():
            sample_id = str(row[args.sample_id_column])
            
            if sample_id in sample_accessions:
                # Make the biosample_accessions the study_accession
                updated_metadata.at[idx, 'sample_accession'] = sample_accessions[sample_id]['biosample_accession']
                updated_metadata.at[idx, 'ena_accession'] = sample_accessions[sample_id]['ena_accession']
        
        # Save the updated metadata
        updated_metadata_path = os.path.join(args.output, 'metadata_with_accessions.tsv')
        updated_metadata.to_csv(updated_metadata_path, sep='\t', index=False)
        logger.info(f"Wrote updated metadata with accessions to {updated_metadata_path}")
        
    else:
        # Create an empty results file
        pd.DataFrame(columns=['sample_alias', 'ena_accession', 'biosample_accession', 'status']) \
          .to_csv(os.path.join(args.output, 'accessions.tsv'), sep='\t', index=False)
        logger.info("No accessions to report")
    
    # Write out a summary report to make things more readable for us
    with open(os.path.join(args.output, 'submission_summary.txt'), 'w') as f:
        f.write("=== ENA XML Submission Summary ===\n")
        f.write(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input metadata file: {args.metadata}\n")
        f.write(f"Study accession: {args.study}\n")
        f.write(f"Environment: {'TEST' if args.test else 'PRODUCTION'}\n\n")
        
        if args.column_mappings and column_mappings:
            f.write(f"Column mappings: Used {len(column_mappings)} mappings from {args.column_mappings}\n\n")
        
        f.write(f"Total samples processed: {len(metadata_df)}\n")
        f.write(f"Successfully submitted: {len(all_results)}\n")
        
        if all_results:
            # Count samples by status
            status_counts = {}
            for result in all_results:
                status = result['status']
                status_counts[status] = status_counts.get(status, 0) + 1
            
            f.write("\nStatus counts:\n")
            for status, count in status_counts.items():
                f.write(f"  {status}: {count}\n")
        
        if all_errors:
            f.write("\nErrors encountered:\n")
            for error in all_errors:
                f.write(f"  {error}\n")
    
    # Return success status for workflow
    success = len(all_results) > 0 and len(all_errors) == 0
    with open(os.path.join(args.output, 'success.txt'), 'w') as f:
        f.write(str(success).lower())
    
    logger.info(f"Submission process complete. Success: {success}")
    
    return 0 if success else 1

if __name__ == "__main__":
    # Will return 0 if successful, 1 if not
    sys.exit(main())