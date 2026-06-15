FROM ubuntu:xenial

# metadata
LABEL base.image="ubuntu:xenial"
LABEL container.version="1"
LABEL software="ARTIC-nCov19"
LABEL software.version="1.1.3"
LABEL description="Conda environment for ARTIC network nCov19 bioinformatic SOP"
LABEL website="https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html"
LABEL license=""
LABEL comments: "Clone of artic ncov2019 at most recent commit as of Jan, 22 2021; Dockerfile based on https://github.com/StaPH-B/docker-builds/blob/master/artic-ncov2019-medaka/Dockerfile"
LABEL maintainer="Kevin Libuit"
LABEL maintainer.email="kevin@libuitsci.com"

# install needed software for conda install
RUN apt-get update && apt-get install -y \
  wget \
  git \
  build-essential

# get miniconda and the artic-ncov19 repo
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
 bash ./Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b  &&\
  rm Miniconda3-latest-Linux-x86_64.sh &&\
  git clone --recursive https://github.com/artic-network/artic-ncov2019.git &&\
  cd artic-ncov2019 &&\
  git checkout b29c6afcdfb7d4ee51425426a9c558c1cbc6b815 

# set the environment
ENV PATH="/miniconda/bin:$PATH"

# create the conda environment and set as default
RUN conda env create -f /artic-ncov2019/environment.yml
RUN echo "source activate artic-ncov2019" > /etc/bash.bashrc
ENV PATH /miniconda/envs/artic-ncov2019/bin:$PATH
WORKDIR /data

