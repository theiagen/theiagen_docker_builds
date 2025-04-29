FROM mambaorg/micromamba:2.0.8 AS app

# build and run as root users since micromamba image has 'mambauser' set as the $USER
USER root
# set workdir to default for building; set to /data at the end
WORKDIR /

LABEL base.image="mambaorg/micromamba:2.0.8"
LABEL dockerfile.version="1"
LABEL description="Various scripts for rsveillance"
LABEL maintainer="Sage Wright"
LABEL maintainer.email="sagemwright@gmail.com"

# Example apt-get command for commonly-missing dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
 wget \
 ca-certificates \
 procps && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Install your desired software into the base conda/micromamba environment, pinning the version
# clean up conda garbage
# make /data to use as a working directory
RUN micromamba install --name base -c conda-forge -c bioconda -c defaults \
 numpy \
 pandas \
 scikit-learn \
 scikit-allel \
 gsl \
 biopython \
 plotnine \
 python=3.12 && \
 micromamba clean -a -y && \
 mkdir /data

# 'ENV' instructions set environment variables that persist from the build into the resulting image
# set the environment, add base conda/micromamba bin directory into path
# set locale settings to UTF-8
ENV PATH="/opt/conda/bin/:${PATH}" \
 LC_ALL=C.UTF-8

WORKDIR /scripts
COPY *.py ./
RUN chmod +x *.py 

ENV PATH="${PATH}:/scripts"

# set final working directory to /data
WORKDIR /data

FROM app AS test

RUN get_depth_distribution.py -h && \
  get_align_stats.py -h && \
  get_depth_windows.py -h && \
  get_pi_ranges.py -h && \
  get_pidp_ranges.py -h && \
  get_refdist_distribution.py -h