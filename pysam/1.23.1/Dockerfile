ARG SAMTOOLS_VER="1.23"

FROM ubuntu:jammy AS builder

ARG SAMTOOLS_VER

# install dependencies required for compiling samtools
RUN apt-get update && apt-get install --no-install-recommends -y \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    zlib1g-dev \
    libssl-dev \
    libdeflate-dev \
    gcc \
    wget \
    make \
    perl \
    bzip2 \
    gnuplot \
    ca-certificates

# download, compile, and install samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER} && \
    ./configure && \
    make && \
    make install && \
    make test

### start of app stage ###
FROM ubuntu:jammy AS app

ARG SAMTOOLS_VER

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="samtools"
LABEL software.version="${SAMTOOLS_VER}"
LABEL description="Tools (written in C using htslib) for manipulating next-generation sequencing data"
LABEL website="https://github.com/samtools/samtools"
LABEL license="https://github.com/samtools/samtools/blob/develop/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

ARG DEBIAN_FRONTEND=noninteractive

# install dependencies required for running samtools
# 'gnuplot' required for plot-ampliconstats
RUN apt-get update && apt-get install --no-install-recommends -y \
    perl \
    zlib1g \
    libncurses5 \
    bzip2 \
    liblzma5 \
    libcurl3-gnutls \
    libdeflate0 \
    gnuplot \
    python3 \
    python3-pip \
    && apt-get autoclean && rm -rf /var/lib/apt/lists/*

# copy in samtools executables from builder stage
COPY --from=builder /usr/local/bin/* /usr/local/bin/

# add local scripts
COPY scripts/*.py /usr/bin/

ENV LC_ALL=C

# final working directory is /data
WORKDIR /data

# install pysam + biopython
RUN pip3 install pysam biopython

# default command is to pull up help options
CMD ["samtools", "--help"]

### start of test stage ###
FROM app AS test

ARG SAMTOOLS_VER

# check PATH
RUN samtools --help

COPY scripts/ /test/scripts/
COPY tests/ /test/tests/

RUN pip3 install pytest

# install make and wget for downloading test files
RUN apt-get update && apt-get install --no-install-recommends -y wget ca-certificates

WORKDIR /test

RUN echo "Testing samtools stats and samtools faidx on SARS-CoV-2 amplicon data" && \
  wget -q https://raw.githubusercontent.com/StaPH-B/docker-builds/master/tests/SARS-CoV-2/SRR13957123.consensus.fa && \
  wget -q https://raw.githubusercontent.com/StaPH-B/docker-builds/master/tests/SARS-CoV-2/SRR13957123.primertrim.sorted.bam && \
  wget -q https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.primer.bed && \
  samtools stats SRR13957123.primertrim.sorted.bam > SRR13957123.stats.txt && \
  echo "Showing first 20 lines of samtools stats output:" && \
  head -n 20 SRR13957123.stats.txt && \
  samtools faidx SRR13957123.consensus.fa && \
  samtools ampliconstats nCoV-2019.primer.bed SRR13957123.primertrim.sorted.bam > SRR13957123_ampliconstats.txt && \
  plot-ampliconstats plot SRR13957123_ampliconstats.txt && \
  python3 /usr/bin/gene_coverage.py \
    --bam SRR13957123.primertrim.sorted.bam \
    --bedfile nCoV-2019.primer.bed && \
  cat COVERAGE_STATS.tsv && \
  ls && \
  pytest -q /test/tests