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

Then to prepare metadat for bulk webin cli:

```
python prepare_ena_data.py \
  --input test_results/metadata_with_accessions.tsv \
  --output ena_prepped_metadata.tsv \ 
  --file-paths file_paths.json \
  --sample-id-column sample_id \
  --study-accession PRJEB12345
```

Then we need to localize the files

```
python localize_files.py --file-path file_paths.json
```

Then test the bulky_webincli.py in the docker:
```
docker run -v ~/phb-dev/theiagen_docker_builds/ena-webinbulk-cli/terra_to_ena:/workdir -it --entrypoint /bin/bash terra_to_ena:0.4
```

And then in the docker:
```
cd /workdir

# Validate first
bulk_webincli.py -u Webin-xxxx  -p  TestPassword -w /scripts/webin-cli.jar -g reads -s ena_prepped_metadata.tsv -d . -m validate -t

# And now submit
bulk_webincli.py -u Webin-xxxx  -p  TestPassword -w /scripts/webin-cli.jar -g reads -s ena_prepped_metadata.tsv -d . -m submit -t
```