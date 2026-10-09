# snpeff container

Main tool: [snpeff](https://pcingola.github.io/SnpEff/)

Full documentation: [SnpEff & SnpSift](https://pcingola.github.io/SnpEff/)

Code repository: [theiagen_docker_builds](https://github.com/theiagen/theiagen_docker_builds)

Derived from the [StaPH-B snpeff 5.4c](https://github.com/StaPH-B/docker-builds/tree/master/build-files/snpeff/5.4c) build.

SnpEff annotates genetic variants and predicts their functional effects (e.g. amino acid
changes) on genes. SnpSift is used after SnpEff annotation to filter and manipulate
annotated files.

Additional tools:
- SnpSift 5.4c (bundled with snpEff)
- snpEff helper scripts, e.g. `buildDbNcbi.sh` (on `PATH` from `/snpEff/scripts`)
- [theiagene](https://github.com/theiagen/theiagene) 1.0.1 (installed via `pip` from the release tarball)

## Example Usage

Print help and confirm the install:

```bash
docker run --rm theiagen/snpeff:5.4c snpeff -version
```

List available databases, then download one and annotate a VCF. Prebuilt databases download
to `/snpEff/data`; mount a host directory there to reuse them across runs:

```bash
docker run --rm theiagen/snpeff:5.4c snpeff databases

mkdir -p $PWD/snpeff_data
docker run --rm \
  -v $PWD:/data \
  -v $PWD/snpeff_data:/snpEff/data \
  theiagen/snpeff:5.4c \
  bash -c "snpeff download -v hg19 && snpeff hg19 /data/input.vcf > /data/annotated.vcf"
```

Build a database from an NCBI GenBank accession (writes `data/` and `snpEff.config` to the
working directory):

```bash
docker run --rm -v $PWD:/data theiagen/snpeff:5.4c buildDbNcbi.sh CP014866.1
```

Filter annotated output with SnpSift:

```bash
docker run --rm -v $PWD:/data theiagen/snpeff:5.4c \
  bash -c 'snpsift filter "(QUAL >= 30)" /data/annotated.vcf > /data/filtered.vcf'
```

Run the bundled `theiagene` CLI:

```bash
docker run --rm theiagen/snpeff:5.4c theiagene --help
docker run --rm theiagen/snpeff:5.4c theiagene report_variants --help
```
