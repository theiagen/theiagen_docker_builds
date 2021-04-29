FROM ubuntu:xenial

LABEL base.image="ubuntu:xenial"
LABEL dockerfile.version="1"
LABEL software="MIDAS"
LABEL software.version="1.0.0"
LABEL description="Rapid genomic-distance comparison for taxonomic identification of microbial pathogens"
LABEL website="https://github.com/hesslab-midas"
LABEL license="NA"
LABEL maintainer="Kevin Libuit"
LABEL maintainer.email="Kevin@Libuit.com"

# install c compiler & needed software for conda install
RUN apt-get update && apt-get install -y \
  wget \
  gcc \
  git \
  build-essential

# get miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
 bash ./Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b  &&\
  rm Miniconda3-latest-Linux-x86_64.sh

# set the environment
ENV PATH="/miniconda/bin:$PATH"\
 LC_ALL=C

# get MIDAS installer and configured_appsdata to image
RUN wget https://storage.googleapis.com/theiagen-public-files/midas_nsphl/configured_appdata.tar.gz &&\
    wget https://storage.googleapis.com/theiagen-public-files/midas_nsphl/midas_installer.tar.gz &&\
    tar -xzf configured_appdata.tar.gz &&\
    tar -xzf midas_installer.tar.gz &&\
    rm configured_appdata.tar.gz &&\
    rm midas_installer.tar.gz
    
# create MIDAS conda environment
RUN conda env create -n "midas" -f /midas_installer/env.yml
SHELL ["conda", "run", "-n", "midas", "/bin/bash", "-c"]
RUN cd /midas_installer/midas && \
  python setup.py build_ext --inplace && \
  python setup.py install && \
  cd ../midas-cli && \
  python setup.py install
ENV PATH /miniconda/envs/midas/bin:$PATH


# initialize MIDAS application
RUN echo "/midas/appdata" | midas init
 

# overwrite appdata directory with configured_appdata
RUN rm -rf /midas/appdata/* && \
  mv /configured_appdata/* /midas/appdata && \
  mkdir /midas/appdata/db/sequences

# set working directory to /data
RUN mkdir /data
WORKDIR /data

