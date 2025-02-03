FROM mambaorg/micromamba:0.27.0 as app

# build and run as root users since micromamba image has 'mambauser' set as the $USER
USER root

# set workdir to default for building; set to /data at the end
WORKDIR /

# ARG variables only persistent for the build stage
ARG MYKROBE_VER="0.12.1"
#Pulling specific commit with Sage's corrections.
ARG SONNEITYPING_VER="7d18a7c"

# metadata labels
LABEL base.image="mambaorg/micromamba:0.27.0"
LABEL dockerfile.version="1"
LABEL software="Mykrobe & Genotyphi & Sonneityping"
LABEL software.version=${MYKROBE_VER}
LABEL description="Conda environment for Mykrobe, particularly for Genotyphi"
LABEL website1="https://github.com/Mykrobe-tools/mykrobe"
LABEL license1="MIT"
LABEL license1.url="https://github.com/Mykrobe-tools/mykrobe/blob/master/LICENSE"
LABEL website2="https://github.com/katholt/genotyphi"
LABEL license2="GNU General Public License v3.0"
LABEL license2.url="https://github.com/katholt/genotyphi/blob/main/LICENSE"
LABEL website3="https://github.com/katholt/sonneityping/"
LABEL license3="GNU General Public License v3.0"
LABEL license3.url="https://github.com/katholt/sonneityping/blob/master/LICENSE.txt"
LABEL maintainer1="Andrew Hale"
LABEL maintainer1.email="andrew.hale@theiagen.com"

# install dependencies; cleanup apt garbage
RUN apt-get update && apt-get install -y --no-install-recommends \
 wget \
 ca-certificates \
 git \
 procps \
 jq && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# get the genotyphi code; make /data
# cloning this commit: 98a6e9ccdf069bb86fcf41035b8c5fa92952aa9e
# url: https://github.com/katholt/genotyphi/commit/98a6e9ccdf069bb86fcf41035b8c5fa92952aa9e
# because genotyphi v1.9.1 does NOT include parse_typhi_mykrobe.py script for parsing mykrobe results
RUN git clone https://github.com/katholt/genotyphi.git && \
 cd genotyphi && \
 git checkout 98a6e9ccdf069bb86fcf41035b8c5fa92952aa9e && \
 chmod +x /genotyphi/parse_typhi_mykrobe.py && \
 mkdir -v /data

# Get the sonneityping code
RUN git clone https://github.com/katholt/sonneityping.git && \
 cd sonneityping && \
 git checkout 7d18a7cbfdb1f931765cc0b243a17605b7cff570 && \
 chmod +x /sonneityping/parse_mykrobe_predict.py

# set the PATH and LC_ALL for singularity compatibility
ENV PATH="${PATH}:/opt/conda/envs/mykrobe/bin/:/genotyphi:/sonneityping" \
 LC_ALL=C.UTF-8

# create the conda environment, install mykrobe via bioconda package; cleanup conda garbage
# INSTALL PANDAS HERE INSTEAD
RUN micromamba create -n mykrobe -y -c conda-forge -c bioconda -c defaults \
 mykrobe=${MYKROBE_VER} \
 python \ 
 pip \
 pandas && \
 micromamba clean -a -y

# so that mamba/conda env is active when running below commands
ENV ENV_NAME="mykrobe"
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# get the latest databases (AKA "panels")
RUN mykrobe panels update_metadata && \
 mykrobe panels update_species all && \
 mykrobe panels describe

WORKDIR /data

# new base for testing
FROM app as test

# so that mamba/conda env is active when running below commands
ENV ENV_NAME="mykrobe"
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# test with TB FASTQs
RUN wget -O test_reads.fq.gz https://ndownloader.figshare.com/files/21059229 && \
 mykrobe predict -s SAMPLE -S tb -o out.json --format json -i test_reads.fq.gz && \
 cat out.json && \
 mykrobe panels describe && \
 mykrobe --version