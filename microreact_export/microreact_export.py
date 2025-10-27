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
  with open(tsv, 'r', newline='', encoding='utf-8') as f:
    df = pd.read_csv(f, sep="\t", dtype=str)

    logger.debug(f"selected_columns: {selected_columns}")
    logger.info(f"Original metadata TSV has {len(df)} records and columns: {df.columns.tolist()}")
    if selected_columns:
      data = df[selected_columns]
    else:
      logger.info("No selected columns provided; using all columns.")
      data = df
    logger.info(f"Parsed metadata TSV with {len(data)} records and columns: {data.columns.tolist()}")
    logger.info(f"First few rows of data:\n{data.head()}")
    
    if date_column:
      logger.info(f"Validating date column: {date_column}")
      df[date_column] = pd.to_datetime(df[date_column], errors='raise')
      logger.info(f"Date column '{date_column}' validated successfully.")

    if remove_file_columns:
      logger.info("Removing columns that contain file URLs.")
      for col in data.columns:
        if data[col].astype(str).str.startswith(('http://', 'https://', 's3://', 'gs://', 'ftp://')).any():
          logger.info(f"Removing column '{col}' as it to contain file URLs.")
          data = data.drop(columns=[col])

    csv_data = data.to_csv(index=False)
    f.close()
  
  headers = csv_data.splitlines()[0].split(',')
  if id_column not in headers:
    logger.error(f"ID column '{id_column}' not found in metadata headers: {headers}")
    logger.info(f"Defaulting to first column '{headers[0]}' as ID column.")
    id_column = headers[0]

  return csv_data, headers, id_column

def encode(data: str) -> str:
  encoded_data = "data:application/octet-stream;base64," + base64.b64encode(data.encode('utf-8')).decode('utf-8')
  file_id = str(uuid.uuid4())
  logger.info(f"File encoded with ID: {file_id}")
  return encoded_data, file_id

def add_trees(tree_files: List[Path], id_column: str):
  print(tree_files)
  tree_files_dict = {}
  tree_dict = {}
  for tree_file in tree_files:
    logger.info(f"Processing tree file: {tree_file}")
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
  return tree_files_dict, tree_dict

def create_project_entry(
  entry_type: str,
  id_column: str,
  remove_file_columns: Optional[bool] = True,
  metadata: Optional[Path] = None, 
  date_column: Optional[str] = None, 
  selected_columns: Optional[List[str]] = None,
  tree_files: Optional[List[Path]] = None
):
  if entry_type == "metadata":
    logger.info("Creating metadata entry")
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
  elif entry_type == "tree":
    logger.info("Creating tree entries")
    tree_files_dict, tree_dict = add_trees(tree_files, id_column)
    return tree_files_dict, tree_dict
  elif entry_type == "map":
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

  metadata_entry, metadata_id, headers = create_project_entry("metadata",  id_column, remove_file_columns, metadata, date_column, selected_columns, None)

  # Create Metadata Dataset and Table Entries
  dataset_id = str(uuid.uuid4())
  logger.info(f"Created dataset entry with ID {dataset_id}, linked to metadata file ID {metadata_id}")
  dataset_entry = {
      "id": dataset_id,
      "file": metadata_id,
      "idFieldName": id_column
  }
  columns = [{"field": col} for col in selected_columns] if selected_columns else [{"field": col} for col in headers]
  table_id = str(uuid.uuid4())
  logger.info(f"Created table entry with ID {table_id}, linked to metadata file ID {metadata_id}")
  table_entry = {
    "title": "Metadata",
    "file": metadata_id,
    "columns": columns
  }    

  if tree_files:
    tree_files_dict, tree_dict = create_project_entry("tree",id_column, remove_file_columns, None, None, None, tree_files)
  else:
    logger.info("No tree files provided; skipping tree entry creation.")
    tree_files_dict, tree_dict = {}, {}
  
  if "latitude" in headers and "longitude" in headers:
    map_entry = create_project_entry("map", "", None, None, None, None)
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
  logger.info(f"Created Microreact project JSON")

  if access_token:
    url = "https://microreact.org/api/projects/create"
    if restricted_access:
      url += "?access=private"
    headers = {
        "Content-Type": "application/json; charset=UTF-8",
        "Access-Token": access_token
    }
    response = requests.post(url, headers=headers, json=project_input)
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
    metadata: Optional[str],
    id_column: Optional[str],
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
  logger.info(f"Fetched current project data for update")
  if metadata:
    metadata_entry, metadata_id, headers = create_project_entry("metadata", id_column, date_column, remove_file_columns, metadata, None, None, None)

    if "latitude" in headers and "longitude" in headers:
      if updated_project.get("maps") is None:
        map_entry = create_project_entry("map", "", None, None, None, None)
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
    tree_files_dict, tree_dict = create_project_entry("tree",id_column, remove_file_columns, None, None, None, tree_files)
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
  argparser.add_argument("--restricted_access", type=bool, default=True, help="Set project access to private if True")
  argparser.add_argument("--update", type=bool, default=False, help="Update an existing Microreact project if True, pair with --project_url")
  argparser.add_argument("--remove_file_columns", type=bool, default=True, help="Remove columns associated with cloud URLs if True")
  argparser.add_argument("--id_column", type=str, help="Column to use as the unique ID field")
  argparser.add_argument("--date_column", type=str, help="Column name for date usage and validation")
  args = argparser.parse_args()

  if args.update and args.project_url:
    microreact_response, updated_project = update_microreact_project(
        project_url=args.project_url,
        access_token=args.access_token,
        metadata=args.metadata_tsv,
        id_column=args.id_column,
        date_column=args.date_column,
        remove_file_columns=args.remove_file_columns,
        tree_files=args.tree_files
    )
    with open("microreact_response.json", "w", encoding="utf-8") as response_file:
      json.dump(microreact_response.json(), response_file, ensure_ascii=False, indent=2)
    with open("project_input.json", "w", encoding="utf-8") as project_file:
      project_file.write(json.dumps(updated_project))
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
    with open("microreact_response.json", "w", encoding="utf-8") as response_file:
      json.dump(microreact_response, response_file, ensure_ascii=False, indent=2)
    with open("project_input.json", "w", encoding="utf-8") as project_file:
      project_file.write(project_json)
if __name__ == "__main__":
  main()
