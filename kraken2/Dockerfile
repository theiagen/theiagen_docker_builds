FROM ubuntu:jammy AS app

# for easy upgrade later. ARG variables only persist during build time.
ARG K2VER="2.17.1"
ARG BRACKENVER="3.1"
# This is supposed to be v2.1.4 but there was a typo on their repo as "v2.14"
LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="Kraken2"
LABEL software.version="${K2VER}"
LABEL description="Taxonomic sequence classifier"
LABEL website="https://github.com/DerrickWood/kraken2"
LABEL license="https://github.com/DerrickWood/kraken2/blob/master/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

# install dependencies and cleanup apt garbage
RUN apt-get update && apt-get -y --no-install-recommends install \
    wget \
    ca-certificates \
    zlib1g-dev \
    make \
    g++ \
    rsync \
    cpanminus \
    ncbi-blast+ && \
    rm -rf /var/lib/apt/lists/* && apt-get autoclean

# perl module required for kraken2-build
RUN cpanm Getopt::Std

# DL Kraken2, unpack, and install
RUN wget https://github.com/DerrickWood/kraken2/archive/v${K2VER}.tar.gz && \
    tar -xzf v${K2VER}.tar.gz && \
    rm -rf v${K2VER}.tar.gz && \
    cd kraken2-${K2VER} && \
    ./install_kraken2.sh . && \
    mkdir /data /kraken2-db

ENV PATH="/kraken2-${K2VER}/:$PATH" \
    LC_ALL=C

# DL Bracken, unpack, and install
RUN wget https://github.com/jenniferlu717/Bracken/archive/refs/tags/v${BRACKENVER}.tar.gz && \
    tar -xzf v${BRACKENVER}.tar.gz && \
    rm -rf v${BRACKENVER}.tar.gz && \
    cd Bracken-${BRACKENVER} && \
    bash install_bracken.sh

ENV PATH="/Bracken-${BRACKENVER}/:/Bracken-${BRACKENVER}/src/:$PATH"

RUN ln -s /usr/bin/python3 /usr/bin/python

# default command is to pull up help options
CMD [ "kraken2", "--help" ]

WORKDIR /data

### start of test stage ###
FROM app AS test 

# print help and versions
RUN kraken2 --help && \
    kraken2-build --help && \
    kraken2 --version && \
    blastn -version && \
    dustmasker -version && \
    segmasker -version && \
    bracken -v