# PulseNet 2.0 Allele Calling container

**Main tool**: [PulseNet 2.0 AlleleCalling](https://github.com/ncezid-biome/pulsenet2.0-bfx/tree/47644186f2df27e9f01a000d47c451135a75f65d/processes/AlleleCalling)

**Code repository**: <https://github.com/ncezid-biome/pulsenet2.0-bfx/tree/47644186f2df27e9f01a000d47c451135a75f65d/processes/AlleleCalling

v1.0.1 of this directory copies the directory linked above exactly with the following changes:

- the [original](https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/processes/AlleleCalling/tests/files/qc_kb/qc.json) `tests/files/qc_kb/qc.json` has been updated to [include all entries for the "AlleleCalling" section of this JSON](https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/main/tests/files/qc_kb/qc.json) as the original one is incomplete.
- the [LICENSE.txt](https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/processes/AlleleCalling/LICENSE.txt), [CHANGELOG.md](https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/processes/AlleleCalling/CHANGELOG.md), and [CONTRIBUTING.md](https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/processes/AlleleCalling/CONTRIBUTING.md) files have been removed from this directory.
- the [original README.md](https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/processes/AlleleCalling/README.md) has been replaced with this one, which is specific to the Docker image used and not the actual tool's usage.

The Dockerfile has been **unchanged** from the original repository.

**Basic information on how to use this tool**:
- executable: `ngs-run AlleleCalling`
- help: `ngs-run AlleleCalling --help`
- description:  the PulseNet 2.0 AlleleCalling algorithm

**Additional information**:

- Run the executable in the `/app` directory for easiest access to the included metadata JSONs, which are in `/app` subdirectories.
- Blast DBs are not included in this docker image but publicly-available ones can be found and used here: `gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/`

**Full documentation**: [See the README.md here](https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/processes/AlleleCalling/README.md)

## Example Usage

This example shows how to analyze a _Campylobacter jejuni_ sample.

1. Mount the appropriate database to the Docker image. In this example, I will mount `"gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/CAMPY.tar.gz"`
2. Ensure your FASTA is gzipped and that the headers are stripped of extra content (`sed -i '/^>/s/[[:space:]].*$//' assembly.fasta && gzip assembly.fasta`)

```bash
cd /app
# uncompress the blast DB
tar -xzf CAMPY.tar.gz

ngs-run AlleleCalling \
  --sample-id my-sample \
  --publish-dir . \
  --n-threads 1 \
  --assembly assembly.fasta.gz \
  --blast-kb.path tests/files/blast_kb \
  --blast-kb.similarity 70.0 \
  --blast-kb.db CAMPY \
  --blast-kb.loci "CAMPY/loci.tsv" \
  --qc-kb.path tests/files/qc_kb \
  --organism.genus CAMPY \
  --organism.species JEJUNI

# organism.genus and organism.species must match the **controlled language** found in the tests/files/qc_kb/qc.json file
```

## Changelog

- v1.0.0 was a direct copy of the original repository and so this version isn't included in this repository
- v1.0.1 has updated the qc.json file as described above, which has led to its inclusion
