# base image
FROM ubuntu:xenial

# metadata
LABEL base.image="ubuntu:xenial"
LABEL description="Theiagen utility box container to contain: Python3, BaseSpace Sequence Hub CLI, gsutil"
LABEL maintainer="Kevin Libuit"
LABEL maintainer.email="Kevin@libuitsci.com"

# install env essentials 
RUN apt-get update && apt-get install -y \
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
  curl


# install basespace sequence hub CLI & basespace copy plugin
RUN mkdir -p /basespace/bin &&\
  wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O /basespace/bin/bs &&\
  chmod u+x /basespace/bin/bs &&\
  wget "https://api.bintray.com/content/basespace/BaseSpace-Copy-BIN/\$latest/linux/bscp?bt_package=develop" -O /basespace/bin/bs-cp     &&\
  chmod u+x /basespace/bin/bs-cp

# export paths
ENV PATH="/basespace/bin:${PATH}"

# install gsutil
RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-327.0.0-linux-x86_64.tar.gz &&\
    tar -xzf google-cloud-sdk-327.0.0-linux-x86_64.tar.gz &&\
    rm google-cloud-sdk-327.0.0-linux-x86_64.tar.gz &&\
    google-cloud-sdk/install.sh --quiet 
    
# export paths
ENV PATH="/google-cloud-sdk/bin:${PATH}"

# set working dir
RUN mkdir /data
WORKDIR /data
