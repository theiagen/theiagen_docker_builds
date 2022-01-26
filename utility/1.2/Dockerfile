FROM ubuntu:focal

# for easy upgrade later; ARG variables do not persist after the image is built
ARG GOOGLE_CLOUD_SDK_VER="370.0.0"
ARG THEIAGEN_UTILITIES_VER="0.1"

# to prevent tzdata from asking for a region during apt updates
ARG DEBIAN_FRONTEND=noninteractive

# metadata
LABEL base.image="ubuntu:focal"
LABEL description="Theiagen utility box container to contain: Python3, BaseSpace Sequence Hub CLI, gsutil"
LABEL maintainer1="Kevin Libuit"
LABEL maintainer1.email="Kevin@libuitsci.com"
LABEL maintainer2="Curtis Kapsak"
LABEL maintainer2.email="curtis.kapsak@theiagen.com"

# install env essentials; remove apt garbage
RUN apt-get update && apt-get install -y --no-install-recommends \
 wget \
 build-essential  \
 libgl1-mesa-glx \
 libegl1-mesa \
 libxrandr2 \
 libxrandr2 \
 libxss1 \
 libxcursor1 \
 libxcomposite1 \
 libasound2 \
 libxi6 \
 libxtst6 \
 python3 \ 
 python3-pip \
 curl \
 zip \
 unzip \
 jq \
 wkhtmltopdf && \
 apt-get autoclean && \
 rm -rf /var/lib/apt/lists/* 

# install basespace sequence hub CLI & basespace copy plugin
RUN mkdir -p /basespace/bin && \
 wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O /basespace/bin/bs && \
 chmod u+x /basespace/bin/bs 

# install gsutil; make /data for use as working directory
RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GOOGLE_CLOUD_SDK_VER}-linux-x86_64.tar.gz && \
 tar -xzf google-cloud-sdk-${GOOGLE_CLOUD_SDK_VER}-linux-x86_64.tar.gz && \
 rm google-cloud-sdk-${GOOGLE_CLOUD_SDK_VER}-linux-x86_64.tar.gz && \
 google-cloud-sdk/install.sh --quiet && \
 mkdir /data 

# install python dependencies
# required, but part of core python: os, sys
RUN pip3 install pdfkit pandas argparse 

# install Theiagen utilities
RUN mkdir /theiagen_utilities && \
 wget https://github.com/theiagen/utilities/archive/refs/tags/v${THEIAGEN_UTILITIES_VER}.tar.gz -P /theiagen_utilities && \
 tar -xzf /theiagen_utilities/v${THEIAGEN_UTILITIES_VER}.tar.gz && \
 rm -v /theiagen_utilities/v${THEIAGEN_UTILITIES_VER}.tar.gz

# set PATH
ENV PATH="/google-cloud-sdk/bin:\
/basespace/bin:\
/theiagen_utilities/scripts:\
${PATH}"

# set working directory as /data
WORKDIR /data
