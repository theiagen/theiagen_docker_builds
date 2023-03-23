# Dockerfile originally written by Rachel St. Jacques and modified by Curtis Kapsak
# original dockerfile here: https://github.com/StaPH-B/docker-builds/blob/master/centroid/1.0.0/Dockerfile

# base image
FROM ubuntu:focal

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="centroid"
LABEL software.version="v1.0.0 + kapsakcj commit"
LABEL description="a tool for determining the ideal reference genome given a set of fasta files"
LABEL original.website="https://github.com/kapsakcj/centroid"
LABEL original.license="https://github.com/kapsakcj/centroid"
LABEL original.website="https://github.com/stjacqrm/centroid"
LABEL original.license="https://github.com/stjacqrm/centroid"
LABEL maintainer="Curtis Kapsak"
LABEL maintainer.email="curtis.kapsak@theiagen.com"

# install dependencies
RUN apt-get update && apt-get -y --no-install-recommends install \
    unzip \
    tar \
    python3 \
    python3-pip \
    gzip \
    wget \
    ca-certificates \
    git && \
    rm -rf /var/lib/apt/lists/* && apt-get autoclean

# install numpy
RUN python3 -m pip install numpy

# install mash and make /data directory
RUN wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar && \
    tar -xvf mash-Linux64-v2.2.tar && \
    rm -v mash-Linux64-v2.2.tar && \
    mkdir /data

# install the centroid script
# pinning to this specific commit: https://github.com/kapsakcj/centroid/commit/6cda1717f6ca34820d6084c52e1f7f4b9620ec9f
RUN git clone https://github.com/kapsakcj/centroid.git && \
    cd centroid && \
    git checkout 6cda171 && \
    chmod 777 centroid.py

#set paths
ENV PATH="/mash-Linux64-v2.2:\
/centroid:\
$PATH"\
    LC_ALL=C

WORKDIR /data

### commenting out since we likely don't need this mash database ###
# RUN mkdir /db && \
#     cd /db && \
#     wget https://gembox.cbcb.umd.edu/mash/RefSeqSketchesDefaults.msh.gz && \
#     gunzip RefSeqSketchesDefaults.msh.gz

# just make sure mash is installed and on PATH
RUN command -v mash && mash --version