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
python3 register_samples_ena.py \
  --metadata test_inputs/bacterial_pathogen_example.tsv \
  --output test_results \
  --username Webin-xxxx \
  --password TestPassword \
  --study PRJEB12345 \
  --sample-id-column sample_id \
  --sample-type prokaryotic_pathogen \
  --test
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