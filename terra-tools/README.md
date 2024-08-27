# Terra Tools Docker image

This docker image contains the following tools:

- Theiagen's [Utilities](https://github.com/theiagen/utilities)
- google cloud SDK (`gsutil`, `gcloud`, etc.)
- Broad Institute's [terra-tools](https://github.com/broadinstitute/terra-tools) scripts (`export_large_tsv.py`, `import_large_tsv.py`, etc.)
- `jq`
- `git`
- `python3` 3.9.2
- pandas 2.0.3
- pdfkit
- pretty-html-table
- scipy
- csvkit (specifically for this command: `csvcut`)

## Docker image tags

Our organization isn't the best. We started off using the version of Theiagen/Utilities to tag the docker images, but then we switched to using dates so we weren't releasing numerous versions for every single code change. It's a little chaotic, sorry.

BUT since dates are more informative to us we are going to stick with using dates. So the date signifies when the docker image was edited, rebuilt, and pushed to the Theiagen Google artifact registry.

## Changelog

`2024-08-27`: added `csvkit` to dockerfile & docker image to allow use of the CLI tool `csvcut`.

## Build & Deploy

Run these commands to build the docker image locally and push it to the Theiagen Google artifact registry:

```bash
# Build the docker image
docker build --progres=plain --tag us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2024-08-27 terra-tools/2024-08-27/

# Push the docker image to the Theiagen Google artifact registry
docker push us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2024-08-27
```
