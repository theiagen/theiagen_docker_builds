# NanoPlot container

Main tool: [NanoPlot](https://github.com/wdecoster/NanoPlot)

Code repository: [NanoPlot](https://github.com/wdecoster/NanoPlot)

Derived from the [StaPH-B NanoPlot 1.48.0 build](https://github.com/StaPH-B/docker-builds/tree/master/build-files/nanoplot/1.48.0), installed through micromamba instead of pip.

Additional tools:
- [SeqKit](https://github.com/shenwei356/seqkit) 2.14.0 - GC content
- environment.yml

## Example Usage

```bash
# NanoPlot summary statistics and HTML report
NanoPlot --fastq reads.fastq.gz --N50 --loglength --tsv_stats --no_static -t 4 -o nanoplot-out

# GC content; GC(%) is column 18 of `seqkit stats -a -T` in SeqKit 2.14.0
seqkit stats -a -T reads.fastq.gz
```

## Changes
Initialize NanoPlot container with SeqKit for GC content
