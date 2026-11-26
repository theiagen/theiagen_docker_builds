# fastq-dl-enafree container

Main tool: [fastq-dl](https://github.com/rpetit3/fastq-dl) 2.0.4

Code repository: [fastq-dl](https://github.com/rpetit3/fastq-dl)

Base image: `us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0`

## Purpose

fastq-dl 2.0.4 queries the European Nucleotide Archive (ENA, `www.ebi.ac.uk`) for run metadata, including with `--provider sra`. From many Google Cloud egress IP addresses the TLS handshake to ENA is reset. fastq-dl 2.0.4 sets no timeout and catches no exception on that request, so the run stops before any download starts.

This image applies a patch to `fastq_dl.py` so that `--provider sra` queries NCBI SRA only. No ENA request is made. The FASTQ download tools (`prefetch`, `fasterq-dump`) are unchanged.

This is a workaround. No upstream fastq-dl release is modified.

## Changes

The patch is `fastq-dl-2.0.4-enafree.patch` (also installed in the image under `/usr/share/doc/fastq-dl-enafree/`). It changes `fastq_dl/fastq_dl.py` in three places:

1. `get_ena_metadata`: the ENA request has a 30 second timeout, and a `requests` exception returns `[False, [None, <message>]]` instead of ending the program.
2. `get_sra_metadata`: `db.fetch_ena_fastq` is replaced by a function that returns an empty `DataFrame`, so pysradb does not query ENA for FASTQ links.
3. `get_run_info`: `only_provider` is true when `--provider sra` is given, so there is no fallback to ENA.

The Dockerfile checks the MD5 of `fastq_dl.py` before the patch (`0a6b79d9a477bbcaafb8a2b4496b9181`) and after it (`053c8c517dcff38a7ce08b3c46c91697`). The build fails if either differs. The patched file has the same MD5 as the one in `us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl-enafree:2.0.4--pyhdfd78af_0-patch1`, the hand-built image that was run on Terra.

## Behavior differences from fastq-dl 2.0.4

- The patch only acts with `--provider sra`. With the default provider or `--provider ena` the image behaves like the base image.
- With `--provider sra`, the `<prefix>-run-info.tsv` file (default prefix `fastq`) has fewer columns because the ENA metadata is not fetched (57 instead of 196 for SRR9207855).
- With `--provider sra`, there is no fallback to downloading from ENA when the SRA download fails.

## Usage

Terra / PHB `fetch_sra_to_fastq`: set the `docker` input to this image and keep `fastq_dl_opts = --provider sra`. To revert, clear the `docker` input.

```bash
docker run --rm -v "$PWD":/data -w /data \
  us-docker.pkg.dev/general-theiagen/theiagen/fastq-dl-enafree:2.0.4-patch1 \
  fastq-dl --accession SRR9207855 --provider sra --outdir out
```

## Tests

The `test` stage prints the version and runs `test_enafree.py`, which needs no network access. It checks the three changes with mocked ENA and SRA calls.

```bash
docker build --target test -t fastq-dl-enafree:test fastq-dl-enafree/2.0.4-patch1/
```

End-to-end check with ENA blocked (needs access to NCBI). The command finishes with exit code 0 and writes `SRR9207855.fastq.gz`, `SRR9207855_1.fastq.gz` and `SRR9207855_2.fastq.gz` to `out/`. The base image fails the same command with a connection error.

```bash
docker build --target app -t fastq-dl-enafree:2.0.4-patch1 fastq-dl-enafree/2.0.4-patch1/

docker run --rm --add-host www.ebi.ac.uk:127.0.0.1 -v "$PWD":/data -w /data \
  fastq-dl-enafree:2.0.4-patch1 \
  fastq-dl --accession SRR9207855 --provider sra --outdir out
```
