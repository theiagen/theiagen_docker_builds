FROM mambaorg/micromamba:2.9.0 AS app

ARG NANOPLOT_VER="1.48.0"
ARG SEQKIT_VER="2.14.0"

LABEL base.image="mambaorg/micromamba:2.9.0"
LABEL dockerfile.version="1"
LABEL software="nanoplot"
LABEL software.version="${NANOPLOT_VER}-seqkit"
LABEL description="Plotting suite for Oxford Nanopore sequencing data and alignments, with SeqKit for GC content"
LABEL website="https://github.com/wdecoster/NanoPlot"
LABEL license="https://github.com/wdecoster/NanoPlot/blob/master/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

# copy version constraints
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# install NanoPlot and SeqKit; --force-pkgs-dirs drops the extracted package cache that --all keeps
RUN micromamba install -y -f /tmp/environment.yml && \
    micromamba clean -y --all --force-pkgs-dirs && \
    rm /tmp/environment.yml

# setup run environment; LC_ALL for singularity compatibility
ENV PATH="/opt/conda/bin:$PATH" \
    CONDA_PREFIX="/opt/conda" \
    LC_ALL=C.UTF-8

USER root

RUN mkdir /data && \
  chown $MAMBA_USER:$MAMBA_USER /data

USER $MAMBA_USER

WORKDIR /data

CMD ["NanoPlot", "--help"]

### start of test stage ###
FROM app AS test

# print help options and versions
RUN NanoPlot --help && \
    NanoPlot --version && \
    seqkit version

# confirm GC content is reported for a small FASTQ
RUN printf "@r1\nACGTGGCC\n+\nIIIIIIII\n" | gzip > test.fastq.gz && \
    seqkit stats -a -T test.fastq.gz | grep -q "75.00"

# download ONT data for a Salmonella isolate, run NanoPlot on the FASTQ file
# Go here for more info: https://www.ebi.ac.uk/ena/browser/view/SRR19787768?show=reads
RUN python3 -c "import urllib.request; urllib.request.urlretrieve('https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/068/SRR19787768/SRR19787768_1.fastq.gz', 'SRR19787768_1.fastq.gz')" && \
    NanoPlot --fastq SRR19787768_1.fastq.gz --loglength --N50 --tsv_stats --no_static -o /data/nanoplot-test && \
    ls -lh /data/nanoplot-test && \
    cat /data/nanoplot-test/NanoStats.txt && \
    grep -q "number_of_reads" /data/nanoplot-test/NanoStats.txt && \
    seqkit stats -a -T SRR19787768_1.fastq.gz
