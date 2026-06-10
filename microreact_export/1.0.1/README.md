# Microreact_Export

`microreact_export.py` creates, updates, and submits Microreact project files in JSON format.

## Inputs

- project_name: Desired name for Microreact project
- project_url: Project ID for current project, used in conjunction with `--update`
- update: Allows for a currently present Microreact project to be updated with new metadata or tree files, used in conjunction with `--project_url`
- metadata_tsv: Input table to be used as Microreact metadata, command line useage only
- tree_files: List of tree files for upload
- selected_columns: Columns to be included in Microreact metadata
- access_token: Access token available on Microreact account page for API access
- restricted_access: Creates a private Microreact project (default: false)
- remove_file_columns: Removes columns associated with cloud filepaths (default: false)
- id_column: Column set to be the unique identifier for Microreact
- date_column: If available, the column name associated with sample dates for Microreact Timeline usage