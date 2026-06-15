# Microreact Tests

`microreact_export_tests.py` is a set of unit tests designed to test `microreact_export.py` versions. The tests are set up to specifically check for functionality with data transformation, Microreact project assembly, and submitting the Microreact projects. More tests are to be added with future added functionality.

## Usage

Export the path of the `microreact_export.py` script you would like to test under the variable `MICROREACT_EXPORT_PATH` and run the test.

## Tests

- test_logger: This test confirms that logger creation works and the logs are written to disk correctly
- test_parse_metadata_tsv_basic: This test creates and parses a simple metadata TSV with no special parameters passed
- test_parse_metadata_tsv_remove_file_columns_and_date: This test performs the metadata parsing with additional parameters of removing columns with file contents and date detection and validation
- test_parse_metadata_tsv_defaults_id_column_if_missing: Test case where a specified `id_column` is not present in the metadata
- test_encode_returns_prefixed_blob_and_uuid: Tests the proper function of encoding 
- test_add_entry_new_and_update: Tests entry addition and updating an existing project JSON with a new entry using `micro_react.add_entry`
- test_create_tree_entry_valid: Tests the tree creation process
- test_create_tree_entry_unsupported_type: Creates a `tree.txt` using an unsupported extension and tests the exception of no valid tree files
- test_create_metadata_entry_with_and_without_timeline: Singularly validates the functionality of timeline creation within `create_metadata_entry`
- test_create_map_entry: Validates the map creation function
- test_create_matrix_entry: Ensures proper creation of SNP matrix entry while also ensuring the removal of `snp-dists` version tag from the matrix
- test_submit_microreact_project_success: Test push to Microreact with access set to restricted
- test_submit_microreact_project_http_error_raises: Test push with 500 HTTP error
- test_update_microreact_project_minimal: Tests a minimal update process with a mock project
- test_update_microreact_project_with_metadata_tree_and_matrix: Tests an update process with metadata, map, tree, and matrix
- test_create_microreact_project_full: Creates a full Microreact project. Includesa map, tree, matrix, and timeline
- test_create_microreact_project_requires_metadata: Tests functionality with missing metadata
- test_main_create_path_writes_project_file: Asserts that the final project file is written to a `.microreact` file
- test_main_update_path_writes_project_file: Asserts that the final project proceeding an update is written to a `.microreact` file