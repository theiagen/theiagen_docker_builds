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
- [theiagene](https://github.com/theiagen/theiagene) NA (installed via `pip` from the cloned repo)

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

Run the bundled `theiagene` CLI:

```bash
docker run --rm theiagen/ensembl-vep:116.0 theiagene --help
docker run --rm theiagen/ensembl-vep:116.0 theiagene gene_coverage --help
```