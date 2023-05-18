# Dockerfile originally written by Rachel St. Jacques and modified by Curtis Kapsak
# original dockerfile here: https://github.com/StaPH-B/docker-builds/blob/master/centroid/1.0.0/Dockerfile

# base image
FROM ubuntu:focal as app

ARG CENTROID_VER="0.1.0"

# metadata
LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="1"
LABEL software="centroid"
LABEL software.version="v0.1.0"
LABEL description="a tool for determining the ideal reference genome given a set of fasta files"
LABEL website="https://github.com/theiagen/centroid"
LABEL license="https://github.com/theiagen/centroid"
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
# ensure mash files are owned by root user and not the UID/GID of 1081147385:1360859114. So odd...
RUN wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar && \
    tar -xvf mash-Linux64-v2.2.tar && \
    rm -v mash-Linux64-v2.2.tar && \
    chown -R root:root /mash-Linux64-v2.2 && \
    mkdir /data

# install the centroid script
RUN wget https://github.com/theiagen/centroid/archive/refs/tags/v${CENTROID_VER}.tar.gz && \
    tar xzf v${CENTROID_VER}.tar.gz && \
    rm -v v${CENTROID_VER}.tar.gz && \
    mv -v centroid-${CENTROID_VER} /centroid && \
    chmod 777 /centroid/centroid.py

#set paths
ENV PATH="/mash-Linux64-v2.2:\
/centroid:\
${PATH}" \
LC_ALL=C

# set final working directory as /data
WORKDIR /data

# testing stage
FROM app as test

# just make sure mash & centroid are installed and on PATH
RUN command -v mash && \
mash --version && \
centroid.py --version && \
centroid.py --help