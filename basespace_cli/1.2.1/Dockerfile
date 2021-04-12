# base image
FROM ubuntu:xenial

# metadata
LABEL base.image="ubuntu:xenial"
LABEL description="BaseSpace Sequence Hub CLI"
LABEL website="https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview"
LABEL maintainer="Kevin Libuit"
LABEL maintainer.email="Kevin@libuitsci.com"

# install env essentials 
RUN apt-get update && apt-get install -y \
  wget \
  build-essential  

# install basespace sequence hub CLI & basespace copy plugin
RUN mkdir -p /basespace/bin &&\
  wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O /basespace/bin/bs &&\
  chmod u+x /basespace/bin/bs &&\
  wget "https://api.bintray.com/content/basespace/BaseSpace-Copy-BIN/\$latest/linux/bscp?bt_package=develop" -O /basespace/bin/bs-cp     &&\
  chmod u+x /basespace/bin/bs-cp

# export paths
ENV PATH="/basespace/bin:${PATH}"

# set working dir
RUN mkdir /data
WORKDIR /data
