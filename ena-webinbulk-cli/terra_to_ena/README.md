# Examples for running 

To run for registring samples:

```
python3 register_samples_ena.py \
  --metadata <path_to_metadata.tsv> \
  --output <output_directory> \
  --username <ena_username> \
  --password <ena_password> \
  --study <study_accession> \
  --sample-id-column <column_name> \
  --sample-type <prokaryotic_pathogen|virus_pathogen>
```