FROM mambaorg/micromamba:2.9.0 AS app

LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"
LABEL description="FASTQ summary statistics with fastq-scan and GC content with SeqKit"
LABEL version="1.0.1-seqkit"

# copy version constraints
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# install fastq-scan, SeqKit, and jq
RUN micromamba install -y -f /tmp/environment.yml && \
    micromamba clean -y --all && \
    rm /tmp/environment.yml

# setup run environment
ENV PATH="/opt/conda/bin:$PATH" \
    CONDA_PREFIX="/opt/conda"

USER root

RUN mkdir /data && \
  chown $MAMBA_USER:$MAMBA_USER /data

USER $MAMBA_USER

WORKDIR /data

### start of test stage ###
FROM app AS test

# print versions and confirm GC content is reported for a small FASTQ
RUN fastq-scan -v && \
    seqkit version && \
    jq --version && \
    printf "@r1\nACGTGGCC\n+\nIIIIIIII\n" | gzip > test.fastq.gz && \
    zcat test.fastq.gz | fastq-scan | jq -e '.qc_stats.read_total == 1' && \
    seqkit stats -a -T test.fastq.gz | grep -q "75.00"
