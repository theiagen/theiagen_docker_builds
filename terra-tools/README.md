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

## Docker image tags

Our organization isn't the best. We started off using the version of Theiagen/Utilities to tag the docker images, but then we switched to using dates so we weren't releasing numerous versions for every single code change. It's a little chaotic, sorry.

BUT since dates are more informative to us we are going to stick with using dates. So the date signifies when the docker image was edited, rebuilt, and pushed to the Theiagen Google artifact registry.

```bash
# example

```
