FROM mambaorg/micromamba:2.5.0-ubuntu24.04 AS app

ARG SEQSENDER_VERSION=1.5.2
ARG BIOFORKLIFT_VERSION=0.5.4
ARG GCLOUD_VERSION=587.0.0

# build and run as root users since micromamba image has 'mambauser' set as the $USER
USER root
# set workdir to default for building; set to /data at the end
WORKDIR /

LABEL base.image="mambaorg/micromamba:2.5.0-ubuntu24.04"
LABEL dockerfile.version="1"
LABEL software="SeqSender"
LABEL software.version="${SEQSENDER_VERSION}"
LABEL description="Python program that automates the process of generating submission files and batch uploading to NCBI archives"
LABEL website="https://github.com/CDCgov/seqsender"
LABEL license="https://github.com/CDCgov/seqsender/blob/master/LICENSE"

# install dependencies and cleanup apt garbage
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    procps \
    libsqlite3-0 && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# download source code; make /data
RUN wget -q https://github.com/CDCgov/seqsender/archive/refs/tags/v${SEQSENDER_VERSION}.tar.gz && \
  tar -xzf v${SEQSENDER_VERSION}.tar.gz && \
  rm -r v${SEQSENDER_VERSION}.tar.gz && \
  mkdir /data

# set PATH and locale settings for singularity
ENV LC_ALL=C.UTF-8 \
  PATH=/opt/conda/bin/:${PATH}:/seqsender-${SEQSENDER_VERSION}

# seqsender pins python 3.9 and its own package versions, so it gets its own env
# run it with `micromamba run -n seqsender seqsender.py ...`
RUN micromamba create --yes --name seqsender -f "/seqsender-${SEQSENDER_VERSION}/env.yaml" && \
  micromamba clean --all --force-pkgs-dirs --yes

# base env: bioforklift (Terra table download/upload) and gcloud (copying reads into the NCBI transfer bucket)
RUN micromamba install --yes --name base -c conda-forge python=3.12 pip biopython=1.85 google-cloud-sdk=${GCLOUD_VERSION} && \
  pip install --no-cache-dir bioforklift==${BIOFORKLIFT_VERSION} && \
  micromamba clean --all --force-pkgs-dirs --yes

# helper scripts used by the PHB SeqSender workflow; they run in the base env
COPY scripts/ /scripts/
RUN chmod +x /scripts/*.py
ENV PATH=${PATH}:/scripts

# default command is to pull up help options
CMD ["micromamba", "run", "-n", "seqsender", "seqsender.py", "--help"]

# set working directory to /data
WORKDIR /data

# test layer
FROM app AS test

# print version and grep for it; if not found, exit with error
RUN micromamba run -n seqsender seqsender.py --help && \
  micromamba run -n seqsender seqsender.py version | \
  grep "Version: ${SEQSENDER_VERSION}" || exit 1

# generate test data for FLU and COV and check if directories are created; otherwise exit with error
RUN for org in FLU COV; \
  do micromamba run -n seqsender seqsender.py test_data -bsn --organism "${org}" --submission_dir test_dir && \
    test -d "test_dir/${org}_TEST_DATA/" || exit 1; \
  done

# test network connection
RUN micromamba run -n seqsender seqsender.py test_network_connection | \
  grep "No network connection issues detected" || exit 1

# check gcloud, bioforklift and the helper scripts
RUN gcloud --version && \
  bioforklift download --help && \
  python3 -c "import pandas, Bio" && \
  build_and_prep.py --help

RUN micromamba list -n base && micromamba list -n seqsender
