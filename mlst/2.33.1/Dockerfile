FROM ubuntu:jammy AS app

ARG MLST_VER="2.33.1"
ARG ANY2FASTA_VER="0.8.1"
ARG BLAST_VER="2.17.0"

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="mlst"
LABEL software.version="${MLST_VER}"
LABEL description="Scan contig files against PubMLST typing schemes"
LABEL website="https://github.com/tseemann/mlst"
LABEL license="https://github.com/tseemann/mlst/blob/master/LICENSE"
LABEL maintainer="Inês Mendes"
LABEL maintainer.email="ines.mendes@theiagen.com"
LABEL maintainer2="Curtis Kapsak"
LABEL maintainer2.email="kapsakcj@gmail.com"

# install dependencies via apt; cleanup apt garbage
# deps needed specifically for db download scripts: libfile-which-perl, curl, parallel
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    libmoo-perl \
    liblist-moreutils-perl \
    libjson-perl \
    libfile-which-perl \
    bioperl \
    gzip \
    file \
    procps \
    zlib1g \
    libgomp1 && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install blast; soft link all binaries to /usr/local/bin
RUN wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VER}/ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz && \
    tar -xzf ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz && \
    rm -v ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz && \
    ln -vs /ncbi-blast-${BLAST_VER}+/bin/* /usr/local/bin/

# get any2fasta; move binary to /usr/local/bin which is already in $PATH
RUN wget https://github.com/tseemann/any2fasta/archive/refs/tags/v${ANY2FASTA_VER}.tar.gz && \
 tar xzf v${ANY2FASTA_VER}.tar.gz && \
 rm v${ANY2FASTA_VER}.tar.gz && \
 chmod +x any2fasta-${ANY2FASTA_VER}/any2fasta && \
 mv -v any2fasta-${ANY2FASTA_VER}/any2fasta /usr/local/bin

# get mlst
RUN wget https://github.com/tseemann/mlst/archive/v${MLST_VER}.tar.gz && \
 tar -xzf v${MLST_VER}.tar.gz && \
 rm -v v${MLST_VER}.tar.gz

# set PATH and perl local settings
ENV PATH="${PATH}:/mlst-${MLST_VER}/bin:" \
    LC_ALL=C.UTF-8

# check dependencies and list available schemes
RUN mlst --check && mlst --list

WORKDIR /data

# default command is to pull up help options; can be overridden of course
CMD ["mlst", "--help"]

### start of test stage ###
# using app stage as base image for test stage
FROM app AS test

# run test
RUN mlst ../mlst-2.33.1/test/example.fna.gz

# need unzip for the datasets tool
RUN apt-get update && apt-get install -y --no-install-recommends unzip

# for pinning the version of NCBI datasets that is downloaded for this test
ARG DATASETS_VER="16.2.0"

# install ncbi datasets tool (pre-compiled binary)
RUN wget -q https://github.com/ncbi/datasets/releases/download/v${DATASETS_VER}/linux-amd64.cli.package.zip && \
  unzip linux-amd64.cli.package.zip && \
  rm linux-amd64.cli.package.zip && \
  chmod +x dataformat datasets

# downloading an E. coli (from CT) with a known newer sequence type ST15199 which was designated on 2023-10-16
# MLST scheme info here (from Enterobase): https://pubmlst.org/bigsdb?page=profileInfo&db=pubmlst_escherichia_seqdef&scheme_id=1&profile_id=15199
# E. coli genome: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032604515.1/
# BioSample: SAMN37796909
ARG GENBANK_ACCESSION="GCA_032604515.1"
RUN datasets download genome accession ${GENBANK_ACCESSION} --filename ${GENBANK_ACCESSION}.zip && \
    mkdir -v ${GENBANK_ACCESSION}-download && \
    unzip ${GENBANK_ACCESSION}.zip -d ${GENBANK_ACCESSION}-download && \
    rm -v ${GENBANK_ACCESSION}.zip && \
    mv -v ${GENBANK_ACCESSION}-download/ncbi_dataset/data/${GENBANK_ACCESSION}/${GENBANK_ACCESSION}*.fna ${GENBANK_ACCESSION}-download/ncbi_dataset/data/${GENBANK_ACCESSION}/${GENBANK_ACCESSION}.genomic.fna && \
    mlst ${GENBANK_ACCESSION}-download/ncbi_dataset/data/${GENBANK_ACCESSION}/${GENBANK_ACCESSION}.genomic.fna | tee mlst.${GENBANK_ACCESSION}.tsv && \
    echo && \
    echo "Checking for ST15199 in the mlst results now..." && \
    grep "15199" mlst.${GENBANK_ACCESSION}.tsv
