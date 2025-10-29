import argparse
import pandas as pd
import requests
import json
import base64
import uuid
from pathlib import Path
from typing import Optional, List
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

def parse_metadata_tsv(tsv: Path, id_column: str, remove_file_columns: bool = True, selected_columns: Optional[List[str]] = None, date_column: Optional[str] = None):
  try:
    if tsv is not None:
      logger.info(f"Parsing metadata TSV file.")
      with open(tsv, 'r', newline='', encoding='utf-8') as f:
        df = pd.read_csv(f, sep="\t", dtype=str)
        logger.info(f"Original metadata TSV has {len(df)} records and columns: {df.columns.tolist()}")
        if selected_columns:
          logger.info(f"Selected_columns: {selected_columns}")
          data = df[selected_columns]
        else:
          logger.info("No selected columns provided; using all columns.")
          data = df
        logger.info(f"Parsed metadata TSV with {len(data)} records and columns: {data.columns.tolist()}")
        logger.info(f"First few rows of data:\n{data.head()}")
        
        if date_column:
          if date_column in df.columns:
            logger.info(f"Validating date column: {date_column}")
            df[date_column] = pd.to_datetime(df[date_column], errors='raise')
            logger.info(f"Date column '{date_column}' validated successfully.")
          else:
            logger.warning(f"Date column '{date_column}' not found in metadata; skipping date validation.")
        else:
          logger.info("No date column provided; skipping date validation.")

        if remove_file_columns:
          logger.info("Removing columns that contain file URLs.")
          for col in data.columns:
            if data[col].astype(str).str.startswith(('http://', 'https://', 's3://', 'gs://', 'ftp://')).any():
              logger.info(f"Removing column '{col}' as it to contain file URLs.")
              data = data.drop(columns=[col])
        else:
          logger.info("Not removing columns containing file URLs.")

        csv_data = data.to_csv(index=False)
        f.close()
  except Exception as e:
    logger.error(f"Error parsing metadata TSV: {e}")
    raise e
  
  headers = csv_data.splitlines()[0].split(',')
  if id_column not in headers:
    logger.error(f"ID column '{id_column}' not found in metadata headers: {headers}")
    logger.info(f"Defaulting to first column '{headers[0]}' as ID column.")
    id_column = headers[0]
  else:
    logger.info(f"ID column found, using specified ID column: {id_column}")

  return csv_data, headers, id_column

def encode(data: str) -> tuple[str, str]:
  encoded_data = "data:application/octet-stream;base64," + base64.b64encode(data.encode('utf-8')).decode('utf-8')
  file_id = str(uuid.uuid4())
  logger.info(f"File encoded with ID: {file_id}")
  return encoded_data, file_id

def add_trees(tree_files: List[Path], id_column: str) -> tuple[dict, dict]:
  tree_files_dict = {}
  tree_dict = {}
  allowed_file_types = {".nwk", ".newick", ".tre", ".tree", ".treefile", ".nhx", ".nex", ".nexus"} 
  for tree_file in tree_files:
    logger.info(f"Processing tree file: {tree_file}")
    if tree_file.suffix.lower() not in allowed_file_types:
      logger.warning(f"Skipping unsupported tree file type: {tree_file.name}")
      continue
    else:
      logger.info(f"Accepted tree file type: {tree_file.name}")
      try:
        with open(tree_file, 'r', encoding='utf-8') as tf:
          tree_content = tf.read()
          tree_encoded, tree_id = encode(tree_content)
          logger.info(f"Tree file encoded with ID: {tree_id}")
          tree_files_dict[tree_id] = {
              "blob": tree_encoded,
              "format": "text/x-nh",
              "file": tree_id,
              "name": tree_file.name,
              "type": "tree"
          }
          tree_dict[f"tree_{tree_files.index(tree_file)+1}"] = {
              "type": "rc",
              "labelField": id_column,
              "title": f"Tree {tree_files.index(tree_file)+1}",
              "showLabels": True,
              "showLeafLabels": True,
              "file": tree_id
          }
      except Exception as e:
        logger.error(f"Error processing tree file {tree_file.name}: {e}")
  return tree_files_dict, tree_dict

def create_metadata_entry(
  id_column: str,
  remove_file_columns: bool = True,
  metadata: Optional[Path] = None, 
  date_column: Optional[str] = None, 
  selected_columns: Optional[List[str]] = None
) -> tuple[dict, str, Optional[List[str]]]:
  logger.info("Creating metadata entry")
  if metadata is None:
    raise ValueError("Metadata file path must be provided for metadata entry creation.")
  metadata_parsed, headers, id_column = parse_metadata_tsv(metadata, id_column, remove_file_columns, selected_columns, date_column)
  metadata_encoded, entry_id = encode(metadata_parsed)
  metadata_entry = {}
  metadata_entry[entry_id] = {
      "blob": metadata_encoded,
      "format": "text/csv",
      "id": entry_id,
      "name": "metadata.csv",
      "type": "data"
  }
  return metadata_entry, entry_id, headers

def create_map_entry() -> dict:
  logger.info("Creating map entry")
  map_entry = {}
  map_entry_id = str(uuid.uuid4())
  map_entry[map_entry_id] = {
      "dataType": "geographic-coordinates",
      "title": "Map",
      "latitudeField": "latitude",
      "longitudeField": "longitude",
      "type": "mapbox",
      "coordinate-unit": "decimal-degrees",
      "showRegions": True,
      "showRegionOutlines": True,
      "regionsColourMethod": "entries",
      "regionsColourOpacity": 100,
      "regionsColourPalette": "ColorBrewer YlOrBr-2"
  }
  return map_entry

def create_microreact_project(
  metadata: Path,
  id_column: str,
  date_column: Optional[str],
  tree_files: Optional[List[Path]],
  access_token: Optional[str],
  restricted_access: bool = True,
  remove_file_columns: bool = True,
  project_name: str = "New Microreact Project",
  selected_columns: Optional[List[str]] = None
):

  metadata_entry, metadata_id, headers = create_metadata_entry(id_column, remove_file_columns, metadata, date_column, selected_columns)

  # Create Metadata Dataset and Table Entries
  dataset_id = str(uuid.uuid4())
  logger.info(f"Created dataset entry with ID {dataset_id}, linked to metadata file ID {metadata_id}")
  dataset_entry = {
    "id": dataset_id,
    "file": metadata_id,
    "idFieldName": id_column
  }
  if selected_columns is not None and headers is not None:
    columns = []
    for col in selected_columns:
      if col in headers:
        logger.info(f"Selected column '{col}' found in metadata headers.")
        columns.append({"field": col})
      else:
        logger.warning(f"Selected column '{col}' not found in metadata headers.")
  else:
    logger.info("No selected columns provided; including all metadata columns in table.")
    columns = [{"field": col} for col in headers] if headers is not None else []

  table_id = str(uuid.uuid4())
  logger.info(f"Created table entry with ID {table_id}, linked to metadata file ID {metadata_id}")
  table_entry = {
    "title": "Metadata",
    "file": metadata_id,
    "columns": columns
  }    

  if tree_files:
    tree_files_dict, tree_dict = add_trees(tree_files, id_column)
  else:
    logger.info("No tree files provided; skipping tree entry creation.")
    tree_files_dict, tree_dict = {}, {}

  if headers is not None:
    if any(field.lower() == "latitude" or field.lower() == "longitude" for field in headers):
      map_entry = create_map_entry()
      logger.info("Map entry created based on presence of latitude and longitude in metadata.")
    else:
      logger.info("No latitude/longitude fields found; skipping map entry creation.")
      map_entry = {}

  project_input = {
    "schema": "https://microreact.org/schema/v1.json",
    "meta": {
      "name": project_name
    },
    "datasets": {
      dataset_id: dataset_entry
    },
    "files": {
      **metadata_entry,
      **tree_files_dict
    },
    "tables": {
      table_id: table_entry
    },
    "trees": {
      **tree_dict
    },
    "charts": {},
    "filters": {},
    "maps": {
      **map_entry
    },
    "networks": {},
    "notes": {},
    "panes": {},
    "slicers": {},
    "styles": {},
    "timelines": {},
    "views": {},
    "matrices": {}
  }

  # Write Input to JSON
  project_input_json = json.dumps(project_input)
  logger.info("Created Microreact project JSON")

  if access_token:
    url = "https://microreact.org/api/projects/create"
    if restricted_access:
      url += "?access=private"
    post_headers = {
      "Content-Type": "application/json; charset=UTF-8",
      "Access-Token": access_token
    }
    response = requests.post(url, headers=post_headers, json=project_input)
    logger.info(f"Microreact API response status: {response.status_code}")
    if response.status_code == 200:
      microreact_response = response.json()
      logger.info(f"Microreact project created successfully with ID: {microreact_response.get('id')}")
    else:
      logger.error(f"Failed to create Microreact project: {response.text}")

  return microreact_response, project_input_json

def update_microreact_project(
  project_url: str,
  access_token: str,
  id_column: str,
  metadata: Optional[Path],
  date_column: Optional[str],
  remove_file_columns: bool,
  tree_files: Optional[List[Path]]
):
  logger.info(f"Updating existing Microreact project at {project_url}")
  url = f"https://microreact.org/api/projects/update?project={project_url}"
  get_url = f"https://microreact.org/api/projects/json?project={project_url}"
  get_request_headers = {
    "Access-Token": access_token
  }
  get_response = requests.get(get_url, headers=get_request_headers)
  updated_project = get_response.json()
  logger.info("Fetched current project data for update")
  if metadata:
    metadata_entry, metadata_id, headers = create_metadata_entry(id_column, remove_file_columns, metadata, date_column, None)

    if any(field.lower() == "latitude" or field.lower() == "longitude" for field in headers):
      if updated_project.get("maps") is None:
        map_entry = create_map_entry()
        updated_project["maps"] = {
          **updated_project.get("maps", {}),
          **map_entry
        }
        logger.info("Updated project with new map entry based on metadata")

    new_dataset_id = str(uuid.uuid4())
    logger.info(f"Created new dataset entry with ID {new_dataset_id}, linked to metadata file ID {metadata_id}")
    updated_project["datasets"] = {
      new_dataset_id: {
        "id": new_dataset_id,
        "file": metadata_id,
        "idFieldName": id_column
      }
    }
    columns = [{"field": col} for col in headers]
    new_table_id = str(uuid.uuid4())
    logger.info(f"Created new table entry with ID {new_table_id}, linked to metadata file ID {metadata_id}")
    updated_project["tables"] = {
      new_table_id: {
        "title": "Metadata",
        "file": metadata_id,
        "columns": columns
      }
    }
    updated_project["files"].update(metadata_entry)
    logger.info("Updated project with new metadata")
  
  if tree_files:
    tree_files_dict, tree_dict = add_trees(tree_files, id_column)
    updated_project["files"].update(tree_files_dict)
    if "trees" not in updated_project:
      updated_project["trees"] = tree_dict
    else:
      updated_project["trees"].update(tree_dict)
    logger.info("Updated project with new tree files")

  post_headers = {
    "Content-Type": "application/json; charset=UTF-8",
    "Access-Token": access_token
  }
  logger.info("Sending updated project data to Microreact API")
  post_response = requests.post(url, headers=post_headers, json=updated_project)
  return post_response, updated_project

def main():
  argparser = argparse.ArgumentParser(description="Create a Microreact project from metadata and tree files.")
  argparser.add_argument("--project_name", type=str, help="Name of the Microreact project")
  argparser.add_argument("--project_url", type=str, help="URL of the Microreact project, used for updating existing projects")
  argparser.add_argument("--metadata_tsv", type=Path, help="Path to the metadata file")
  argparser.add_argument("--tree_files", nargs="*", type=Path, help="Paths to the tree files")
  argparser.add_argument("--selected_columns", nargs="*", type=str, help="Columns to include in the Microreact table")
  argparser.add_argument("--access_token", type=str, help="Access token for Microreact API")
  argparser.add_argument("--restricted_access", action="store_true", help="Set project access to private")
  argparser.add_argument("--update", action="store_true", help="Update an existing Microreact project, pair with --project_url")
  argparser.add_argument("--remove_file_columns", action="store_true", help="Remove columns associated with cloud URLs")
  argparser.add_argument("--id_column", type=str, help="Column to use as the unique ID field")
  argparser.add_argument("--date_column", type=str, help="Column name for date usage and validation")
  args = argparser.parse_args()

  if args.update and args.project_url:
    microreact_response, project_json = update_microreact_project(
        project_url=args.project_url,
        access_token=args.access_token,
        metadata=args.metadata_tsv,
        id_column=args.id_column,
        date_column=args.date_column,
        remove_file_columns=args.remove_file_columns,
        tree_files=args.tree_files
    )
  else:
    microreact_response, project_json = create_microreact_project(
      metadata=args.metadata_tsv,
      tree_files=args.tree_files,
      access_token=args.access_token,
      restricted_access=args.restricted_access,
      project_name=args.project_name,
      id_column=args.id_column,
      date_column=args.date_column,
      remove_file_columns=args.remove_file_columns,
      selected_columns=args.selected_columns
    )
  try:
    if microreact_response is not None:
      if hasattr(microreact_response, "json"):
        microreact_response_json = microreact_response.json()
      else:
        microreact_response_json = microreact_response
      with open("microreact_response.json", "w", encoding="utf-8") as response_file:
        json.dump(microreact_response_json, response_file, ensure_ascii=False, indent=2)
    else:
      logger.error("Microreact response is None; cannot write to file.")
  except Exception as e:
    logger.error(f"Error writing Microreact response to file: {e}")

  try:
    if project_json is not None:
      with open("project_input.json", "w", encoding="utf-8") as project_file:
        project_file.write(json.dumps(project_json))
    else:
      logger.error("Updated project is None; cannot write to file.")
  except Exception as e:
    logger.error(f"Error writing updated project to file: {e}")

if __name__ == "__main__":
  main()
