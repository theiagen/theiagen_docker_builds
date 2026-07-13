# ensembl-vep container

Main tool: [ensembl-vep](https://github.com/Ensembl/ensembl-vep)

Full documentation: [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)

Code repository: [theiagen_docker_builds](https://github.com/theiagen/theiagen_docker_builds)

VEP (Variant Effect Predictor) determines the effect of your variants (SNPs, insertions,
deletions, CNVs, or structural variants) on genes, transcripts, protein sequence, and
regulatory regions.

Additional tools:
- htslib / tabix / bgzip (installed via `INSTALL.pl`)
- filter_vep, haplo, variant_recoder (bundled with ensembl-vep)

## What this image is

This is a **base** image: it installs the ensembl-vep API plus Faidx/htslib pinned to the
`release/116.0` git tag (`perl INSTALL.pl --AUTO a`). To keep the image small and generic, it
does **not** bundle any annotation data. The following are intentionally left out and should be
downloaded or mounted at runtime:

- **Cache files** (`--AUTO c`) — the recommended annotation source for offline use
- **FASTA files** (`--AUTO f`)
- **Plugins** (`--AUTO p` / `--PLUGINS`)

VEP itself lives in `/opt/vep/src/ensembl-vep` and is on the `PATH`.

## Building the image

```bash
# Build the production image exactly as CI does (from repo root)
docker build --target app -t theiagen/ensembl-vep:116.0 ./ensembl-vep/116.0

# Build and run the test stage (tests run as RUN steps; a successful build == passing tests)
docker build --target test -t ensembl-vep:test ./ensembl-vep/116.0
```

## Example Usage

Print help and confirm the install:

```bash
docker run --rm theiagen/ensembl-vep:116.0 vep --help
```

Run VEP against a cache. Download a cache once, then mount it read-only at `/opt/vep/.vep`
(VEP's default cache location is `$HOME/.vep`, and `$HOME` is `/root` in this image, so point
`--dir_cache` at the mounted path):

```bash
# one-time: download and unpack a species/assembly cache on the host
mkdir -p $PWD/vep_cache
cd $PWD/vep_cache
curl -O https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz
tar -xzf homo_sapiens_vep_116_GRCh38.tar.gz
cd -

# annotate a VCF fully offline using the mounted cache
docker run --rm \
  -v $PWD:/data \
  -v $PWD/vep_cache:/opt/vep/.vep \
  theiagen/ensembl-vep:116.0 \
  vep --input_file /data/input.vcf \
      --output_file /data/output.vcf \
      --vcf \
      --cache --offline \
      --dir_cache /opt/vep/.vep \
      --species homo_sapiens --assembly GRCh38 \
      --force_overwrite
```

Filter VEP output:

```bash
docker run --rm -v $PWD:/data theiagen/ensembl-vep:116.0 \
  filter_vep --input_file /data/output.vcf --filter "IMPACT is HIGH"
```

## Notes

- Base image: `ubuntu:jammy` (Perl 5.34, which satisfies VEP's Perl >= 5.22 requirement).
- Pinned to the upstream `release/116.0` git tag for reproducibility.
- Bundled Perl dependencies include `DBD::mysql` (database/cache access), `Set::IntervalTree`
  (Haplosaurus + performance), `JSON` (`--json` output), and `PerlIO::gzip` (gzipped input).
- To also install plugins, cache, or FASTA into a derived image, re-run `INSTALL.pl` with the
  relevant `--AUTO` letters (`c`, `f`, `p`) — see the
  [installation docs](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html).
