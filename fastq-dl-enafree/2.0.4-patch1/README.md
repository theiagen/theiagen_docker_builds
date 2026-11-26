# fastq-dl-enafree container

Main tool: [fastq-dl](https://github.com/rpetit3/fastq-dl) 2.0.4

Code repository: [fastq-dl](https://github.com/rpetit3/fastq-dl)

Base image: `us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0`

## Purpose

Hot fix for the PHB `fetch_sra_to_fastq` task (`fastq_dl_sra`), which fails when the connection to the European Nucleotide Archive (ENA, `www.ebi.ac.uk`) is reset.

- PHB runs `fastq-dl` with `--provider sra` by default. The FASTQ files come from NCBI SRA (`prefetch`, `fasterq-dump`), not from ENA.
- fastq-dl 2.0.4 still contacts ENA for run metadata, including with `--provider sra`. From many Google Cloud egress IP addresses the TLS handshake to ENA is reset. fastq-dl 2.0.4 sets no timeout and catches no exception on that request, so the run stops before any download starts. The error shown by the task is `zcat: <accession>_1.fastq.gz: No such file`.
- This image handles a failed ENA request without stopping, and with `--provider sra` it makes no ENA request at all. The FASTQ download from NCBI is the same as before.

The hot fix is limited to the container:

- No WDL change. The PHB task and workflow code are not modified.
- It is applied through the Terra `docker` input of `fetch_sra_to_fastq`. Clearing the input reverts to the original image.
- `fastq_dl_opts` stays `--provider sra`.
- No upstream fastq-dl release is modified.

## Changes

The patch is `fastq-dl-2.0.4-enafree.patch` (also installed in the image under `/usr/share/doc/fastq-dl-enafree/`). It changes `fastq_dl/fastq_dl.py` in three places:

| # | Function | Change | Effect |
| --- | --- | --- | --- |
| 1 | `get_ena_metadata` | The ENA request has a 30 second timeout and is wrapped in `try/except requests.exceptions.RequestException`. The handler returns `[False, [None, <message>]]`. | A failed ENA request is reported as an unsuccessful query instead of ending the program. |
| 2 | `get_sra_metadata` | `db.fetch_ena_fastq` is replaced by a function that returns an empty `DataFrame`. | pysradb does not send its own ENA request for FASTQ links. |
| 3 | `get_run_info` | `only_provider` is true when `--provider sra` is given. | With `--provider sra` there is no ENA metadata query and no fallback to ENA. |

The try/except in change 1 covers the request that fastq-dl makes itself. The ENA request made inside pysradb has no exception handler, so with `--provider sra` changes 2 and 3 are what prevent it.

The Dockerfile checks the MD5 of `fastq_dl.py` before the patch (`0a6b79d9a477bbcaafb8a2b4496b9181`) and after it (`053c8c517dcff38a7ce08b3c46c91697`). The build fails if either differs. The patched file has the same MD5 as the one in `us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl-enafree:2.0.4--pyhdfd78af_0-patch1`, the hand-built image that was run on Terra.

## Impact of removing ENA

With `--provider sra` (the PHB default):

- FASTQ download is unchanged: same tools, same NCBI SRA source, same output file names.
- No ENA request is made, so an ENA connection reset no longer affects the task.
- The `<prefix>-run-info.tsv` file (default prefix `fastq`, returned by PHB as `fastq_metadata`) has fewer columns because the ENA metadata is not fetched: 57 instead of 196 for SRR9207855. PHB does not read this file. Anything that reads `fastq_metadata` downstream sees the smaller file.
- There is no fallback to downloading from ENA when the SRA download fails.

With the default provider or `--provider ena`:

- The ENA request can still be made. It has a 30 second timeout and its failure no longer ends the program.
- Local check with ENA blocked and `fastq-dl --accession SRR9207855 --only-download-metadata -m 2` (no `--provider`, the same shape as the PHB `fetch_srr_accession` command): the base image exits with code 1 and a `ConnectionError`. This image waits 10 seconds after each failed ENA attempt, then queries NCBI SRA, exits with code 0 and writes a 57-column `fastq-run-info.tsv`.
- Not tested: FASTQ download with these providers, and `--provider ena`.

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

Metadata-only check with the default provider and ENA blocked. The command exits with code 0 and writes `out/fastq-run-info.tsv`. The base image exits with code 1.

```bash
docker run --rm --add-host www.ebi.ac.uk:127.0.0.1 -v "$PWD":/data -w /data \
  fastq-dl-enafree:2.0.4-patch1 \
  fastq-dl --accession SRR9207855 --only-download-metadata -m 2 --outdir out
```
