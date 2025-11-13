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

json_scheme = {
  "schema": "https://microreact.org/schema/v1.json",
  "meta": {
    "name": {}
  },
  "datasets": {},
  "files": {},
  "tables": {},
  "trees": {},
  "charts": {},
  "filters": {},
  "matrices": {},
  "maps": {},
  "networks": {},
  "notes": {},
  "panes": {},
  "slicers": {},
  "styles": {},
  "timelines": {},
  "views": {}
}

def parse_metadata_tsv(tsv: Path, id_column: str, remove_file_columns: bool = True, selected_columns: Optional[List[str]] = None, date_column: Optional[str] = None):
  try:
    if tsv is not None:
      dates_validated = False
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
            df[f"{date_column}_validated"] = pd.to_datetime(df[date_column], errors='coerce')
            logger.info(f"Date column '{date_column}' validated successfully.")
            dates_validated = True
            if df[f"{date_column}_validated"].isnull().any():
              logger.warning(f"Some entries in date column '{date_column}' could not be parsed and were set to NaT in {date_column}_validated.")
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

  return csv_data, headers, id_column, dates_validated

def encode(data: str) -> tuple[str, str]:
  encoded_data = "data:application/octet-stream;base64," + base64.b64encode(data.encode('utf-8')).decode('utf-8')
  file_id = str(uuid.uuid4())
  logger.info(f"File encoded with ID: {file_id}")
  return encoded_data, file_id

def add_entry(project: dict, section: str, entry_data: dict) -> dict:
  if section not in project:
    project[section] = entry_data
    logger.info(f"Added new section '{section}' to project.")
  else:
    project[section].update(entry_data)
    logger.info(f"Updated existing section '{section}' in project.")
  return project

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
  metadata_parsed: str, 
  dates_validated: str,
  date_column: Optional[str]
) -> tuple[dict, dict, str]:
  logger.info("Creating metadata entry")
  metadata_encoded, entry_id = encode(metadata_parsed)
  metadata_entry = {}
  metadata_entry[entry_id] = {
      "blob": metadata_encoded,
      "format": "text/csv",
      "id": entry_id,
      "name": "metadata.csv",
      "type": "data"
  }

  timeline_dict = {}
  if dates_validated:
    timeline_dict[date_column] = {
      "controls": False,
      "nodeSize": 14,
      "playing": False,
      "speed": 1,
      "style": "bar",
      "title": date_column,
      "paneId": date_column,
      "dataType": "formatted-value",
      "valueField": f"{date_column}_validated"
    }

  return metadata_entry, timeline_dict, entry_id

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

def create_matrix_entry(
  matrix_file: Path
) -> tuple[dict, dict]:
  try:
    if matrix_file is not None:
      with open(matrix_file, 'r', encoding='utf-8') as mf:
        matrix_content = mf.read()
        matrix_encoded, matrix_id = encode(matrix_content)
        logger.info(f"Matrix file encoded with ID: {matrix_id}")
        matrix_file_entry = {}
        matrix_entry = {}
        logger.info("Creating matrix file entry")
        matrix_file_entry[matrix_id] = {
            "blob": matrix_encoded,
            "format": "text/csv",
            "id": matrix_id,
            "name": matrix_file.name,
            "type": "matrix"
        }
        logger.info("Creating matrix entry")
        matrix_entry[f"matrix_{matrix_id}"] = {
          "controls": True,
          "labelsFontSize": 12,
          "axisLabelsFontSize": 12,
          "showLabels": False,
          "truncateLabels": True,
          "rotateAxisLabels": 90,
          "title": "Matrix",
          "paneId": "matrix-panel",
          "file": matrix_id
        }
        logger.info("Matrix entries created successfully")
  except Exception as e:
    logger.error(f"Error creating matrix entry: {e}")
    raise e
  return matrix_file_entry, matrix_entry

def submit_microreact_project(
  project_input: dict,
  access_token: str,
  restricted_access: bool = True
) -> dict:
  if access_token:
    logger.info("Access token provided, submitting Microreact project via API")
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

  return microreact_response

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
    metadata_parsed, headers, id_column, dates_validated = parse_metadata_tsv(metadata, id_column, remove_file_columns, None, date_column)
    metadata_entry, timeline_dict, metadata_id = create_metadata_entry(metadata_parsed, dates_validated, date_column)

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

  if date_column:
    updated_project["timelines"].update(timeline_dict)

  post_headers = {
    "Content-Type": "application/json; charset=UTF-8",
    "Access-Token": access_token
  }
  logger.info("Sending updated project data to Microreact API")
  post_response = requests.post(url, headers=post_headers, json=updated_project)
  return post_response, json.dumps(updated_project)

def create_microreact_project(
  metadata: Path,
  id_column: str,
  matrix_file: Optional[Path],
  date_column: Optional[str],
  tree_files: Optional[List[Path]],
  remove_file_columns: bool = True,
  project_name: str = "New Microreact Project",
  selected_columns: Optional[List[str]] = None,
  json_scheme: dict = json_scheme
):
  json_scheme["meta"]["name"] = project_name
  logger.info(f"Set project name to: {project_name}")
  if metadata is None:
    raise ValueError("Metadata file path must be provided for metadata entry creation.")
  metadata_parsed, headers, id_column, dates_validated = parse_metadata_tsv(metadata, id_column, remove_file_columns, selected_columns, date_column)
  metadata_entry, timeline_dict, metadata_id = create_metadata_entry(metadata_parsed, dates_validated, date_column)
  json_scheme = add_entry(json_scheme, "files", metadata_entry)
  json_scheme = add_entry(json_scheme, "timelines", timeline_dict)
  
  # Create Metadata Dataset and Table Entries
  dataset_id = str(uuid.uuid4())
  logger.info(f"Created dataset entry with ID {dataset_id}, linked to metadata file ID {metadata_id}")
  dataset_entry = {
    "id": dataset_id,
    "file": metadata_id,
    "idFieldName": id_column
  }
  json_scheme = add_entry(json_scheme, "datasets", {dataset_id: dataset_entry})
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
  json_scheme = add_entry(json_scheme, "tables", {table_id: table_entry})
  if tree_files:
    tree_files_dict, tree_dict = add_trees(tree_files, id_column)
    json_scheme = add_entry(json_scheme, "files", tree_files_dict)
    json_scheme = add_entry(json_scheme, "trees", tree_dict)
  else:
    logger.info("No tree files provided; skipping tree entry creation.")

  if matrix_file is not None:
    logger.info(f"Matrix file provided: {matrix_file.name}")
    matrix_file_entry, matrix_entry = create_matrix_entry(matrix_file)
    json_scheme = add_entry(json_scheme, "files", matrix_file_entry)
    json_scheme = add_entry(json_scheme, "matrices", matrix_entry)
  else:
    logger.info("No matrix file provided; skipping matrix entry creation.")

  if headers is not None:
    if any("latitude" in field.lower() or "longitude" in field.lower() for field in headers):
      map_entry = create_map_entry()
      json_scheme = add_entry(json_scheme, "maps", map_entry)
      logger.info("Map entry created based on presence of latitude and longitude in metadata.")
    else:
      logger.info("No latitude/longitude fields found; skipping map entry creation.")

  # Write Input to JSON
  project_input_json = json.dumps(json_scheme, ensure_ascii=False, indent=2)
  logger.info("Created Microreact project JSON")

  return project_input_json

def main():
  argparser = argparse.ArgumentParser(description="Create a Microreact project from metadata and tree files.")
  argparser.add_argument("--project_name", type=str, help="Name of the Microreact project")
  argparser.add_argument("--project_url", type=str, help="URL of the Microreact project, used for updating existing projects")
  argparser.add_argument("--metadata_tsv", type=Path, help="Path to the metadata file")
  argparser.add_argument("--matrix_file", type=Path, help="Path to the distance matrix file")
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
    project_json = create_microreact_project(
      metadata=args.metadata_tsv,
      matrix_file=args.matrix_file,
      tree_files=args.tree_files,
      project_name=args.project_name,
      id_column=args.id_column,
      date_column=args.date_column,
      remove_file_columns=args.remove_file_columns,
      selected_columns=args.selected_columns,
      json_scheme=json_scheme
    )
    if args.access_token:
      microreact_response = submit_microreact_project(
        project_input=json.loads(project_json),
        access_token=args.access_token,
        restricted_access=args.restricted_access
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
    else:
      logger.info("No access token provided; skipping project submission.")   
  
  try:
    if project_json is not None:
      with open(f"{args.project_name}_input.microreact", "w", encoding="utf-8") as project_file:
        project_file.write(project_json)
    else:
      logger.error("Updated project is None; cannot write to file.")
  except Exception as e:
    logger.error(f"Error writing updated project to file: {e}")

if __name__ == "__main__":
  main()
