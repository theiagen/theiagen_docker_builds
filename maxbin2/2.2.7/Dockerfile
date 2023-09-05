FROM ubuntu:xenial

LABEL base.image="ubuntu:xenial"
LABEL dockerfile.version="1"
LABEL software="MaxBin2"
LABEL software.version="2.2.7"
LABEL description="An automatic tool for binning metagenomics sequences"
LABEL website="https://sourceforge.net/projects/maxbin2"
LABEL license="NA"
LABEL maintainer="Ines Mendes"
LABEL maintainer.email="ines.mendes@theiagen.com"

# install c compiler & needed software for source install
RUN apt-get update && apt-get install -y \
  wget \
  gcc \
  git \
  build-essential

WORKDIR /NGStools

# get miniconda - python2
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN chmod +x Miniconda2-latest-Linux-x86_64.sh
RUN ./Miniconda2-latest-Linux-x86_64.sh -b -p /NGStools/miniconda
ENV PATH="/NGStools/miniconda/bin:$PATH"

# configure conda channels
RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda

# INSTALL maxbin2
RUN conda install maxbin2=2.2.7

# set working directory to /data
WORKDIR /data