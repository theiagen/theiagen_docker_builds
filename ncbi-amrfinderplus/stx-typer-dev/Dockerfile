FROM ubuntu:jammy AS app

ARG AMRFINDER_VER="3.12.8"
ARG AMRFINDER_DB_VER="2024-05-02.2"
ARG BLAST_VER="2.15.0"

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="NCBI AMRFinderPlus"
LABEL software.version="${AMRFINDER_VER}"
LABEL description="NCBI resistance gene detection tool"
LABEL website="https://github.com/ncbi/amr"
LABEL license="https://github.com/ncbi/amr/blob/master/LICENSE"
LABEL maintainer="Curtis Kapsak"
LABEL maintainer.email="kapsakcj@gmail.com"


# ncbi-blast+ installed via apt is v2.12.0 - DISABLING so that we can manually install v2.14.0
# see here for reason why I'm manualy installing 2.14.0 instead of using apt-get: https://github.com/ncbi/amr/releases/tag/amrfinder_v3.11.8

# hmmer installed via apt is v3.3.2
# libgomp1 required for makeblastdb
RUN apt-get update && apt-get install -y --no-install-recommends \
  ca-certificates \
  wget \ 
  curl \
  libgomp1 \
  hmmer \
  procps \
  gzip \
  make \
  gcc \
  g++ \
  git \
  libcurl4 \
  libcurl4-openssl-dev && \
  apt-get autoclean && \
  rm -rf /var/lib/apt/lists/*

# download and install amrfinderplus pre-compiled binaries; make /data
# RUN mkdir amrfinder && cd /amrfinder && \
#   echo "downloading amrfinderplus v${AMRFINDER_VER} pre-compiled binaries from GitHub..." && \
#   wget -q https://github.com/ncbi/amr/releases/download/amrfinder_v${AMRFINDER_VER}/amrfinder_binaries_v${AMRFINDER_VER}.tar.gz && \
#   tar zxf amrfinder_binaries_v${AMRFINDER_VER}.tar.gz && \
#   rm amrfinder_binaries_v${AMRFINDER_VER}.tar.gz && \
#   mkdir /data

# install ncbi-blast linux binaries
RUN echo "downloading ncbi-blast-${BLAST_VER}+ linux binaries from NCBI FTP..." && \
  wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VER}/ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz && \
  tar -xzf ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz && \
  rm -v ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz

# set PATH and locale settings for singularity compatibiliity, set amrfinder and manually-installed blast as higher priority in PATH
ENV PATH="/amrfinder/amr:/ncbi-blast-${BLAST_VER}+/bin:$PATH" \
  LC_ALL=C

# blast has to be installed first so that 'makeblastdb' is available for amrfinder during database update
# download and install amrfinderplus from source; make /data
# TODO - pin to specific commit of ncbi/amr repo?
RUN mkdir amrfinder && cd /amrfinder && \
  git clone https://github.com/ncbi/amr.git && \
  cd amr && \
  git checkout stxtype && \
  git submodule update --init --recursive && \
  make && \
  ./amrfinder -u && \
  make test && \
  mkdir -v /data

# download databases and index them
# done in this manner to pin the database version instead of pulling the latest version with `amrfinder -u` 
# softlink is required for `amrfinder -l` and typical `amrfinder` use cases to work properly
# RUN mkdir -p /amrfinder/data/${AMRFINDER_DB_VER} && \
#   wget -q -P /amrfinder/data/${AMRFINDER_DB_VER} ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/${AMRFINDER_DB_VER}/* && \
#   amrfinder_index /amrfinder/data/${AMRFINDER_DB_VER} && \
#   ln -s /amrfinder/data/${AMRFINDER_DB_VER} /amrfinder/data/latest

# set final working directory
WORKDIR /data

# default command is to print help options
CMD [ "amrfinder", "--help" ]

## Test stage
FROM app AS test

# list database version and available --organism options
RUN amrfinder -l

# run recommended tests from amrfinder
RUN amrfinder --threads 1 --plus -p /amrfinder/amr/test_prot.fa -g  /amrfinder/amr/test_prot.gff -O Escherichia --print_node > test_prot.got && \
  diff /amrfinder/amr/test_prot.expected test_prot.got && \
  amrfinder --threads 1 --plus -n /amrfinder/amr/test_dna.fa -O Escherichia --print_node > test_dna.got && \
  diff /amrfinder/amr/test_dna.expected test_dna.got && \
  amrfinder --threads 1 --plus -n /amrfinder/amr/test_dna.fa -p /amrfinder/amr/test_prot.fa -g /amrfinder/amr/test_prot.gff -O Escherichia --print_node > test_both.got && \
  diff /amrfinder/amr/test_both.expected test_both.got

# run amrfinder on Salmonella, without and with --organism option
RUN wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/941/835/GCA_010941835.1_PDT000052640.3/GCA_010941835.1_PDT000052640.3_genomic.fna.gz  && \
  gzip -d GCA_010941835.1_PDT000052640.3_genomic.fna.gz && \
  amrfinder --threads 1 --plus --nucleotide GCA_010941835.1_PDT000052640.3_genomic.fna --output test1.txt && \
  amrfinder --threads 1 --plus --nucleotide GCA_010941835.1_PDT000052640.3_genomic.fna --organism Salmonella --output test2.txt && \
  cat test1.txt test2.txt

# run amrfinder on Klebesiella oxytoca using --organism/-O flag
RUN wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/812/925/GCA_003812925.1_ASM381292v1/GCA_003812925.1_ASM381292v1_genomic.fna.gz && \
  gzip -d GCA_003812925.1_ASM381292v1_genomic.fna.gz && \
  amrfinder --threads 1 --plus --name GCA_003812925.1 -n GCA_003812925.1_ASM381292v1_genomic.fna  -O Klebsiella_oxytoca -o GCA_003812925.1-amrfinder.tsv 
