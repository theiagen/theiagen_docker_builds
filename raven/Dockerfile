FROM mambaorg/micromamba:2.0.3-ubuntu22.04 AS app

# build and run as root users since micromamba image has 'mambauser' set as the $USER
USER root
# set workdir to default for building; set to /data at the end
WORKDIR /

ARG RAVEN_VER="1.8.3"

# Label the image with metadata that might be important to the user
LABEL base.image="mambaorg/micromamba:2.0.3-ubuntu22.04"
LABEL dockerfile.version="1"
LABEL software="raven"
LABEL software.version="${RAVEN_VER}"
LABEL description="This software runs Raven long read assembler."
LABEL website="https://github.com/theiagen/theiagen_docker_builds/"
LABEL license="https://github.com/theiagen/theiagen_docker_builds/blob/master/LICENSE"
LABEL maintainer="Theron James"
LABEL maintainer.email="theron.james@theiagen.com"

# otherwise tzdata won't install
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install --no-install-recommends -y \
  wget \
  ca-certificates \
  procps && \
  apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN micromamba install --yes --name base -c conda-forge -c bioconda \
  raven-assembler=${RAVEN_VER} && \
  micromamba clean -a -y

ENV PATH="/opt/conda/bin/:${PATH}" \
  LC_ALL=C.UTF-8

# final working directory in "app" layer is /data for passing data in/out of container
WORKDIR /data
