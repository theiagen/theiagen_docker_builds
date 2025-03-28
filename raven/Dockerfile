FROM ubuntu:jammy as app

ARG RAVEN_VER="1.8.3"

# Label the image with metadata that might be important to the user
LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="SoftwareName"
LABEL software.version="${SOFTWARENAME_VER}"
LABEL description="This software runs Raven long read assembler."
LABEL website="https://github.com/theiagen/theiagen_docker_builds/"
LABEL license="https://github.com/theiagen/theiagen_docker_builds/blob/master/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

# Install dependencies for building Raven
RUN apt-get update && apt-get install -y --no-install-recommends \
   wget \
   ca-certificates \
   cmake \
   make \
   git \
   gcc \
   zlib1g-dev \
   build-essential && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Download and extract the Raven release
RUN wget https://github.com/lbcb-sci/raven/archive/refs/tags/${RAVEN_VER}.tar.gz && \
   tar -xzvf ${RAVEN_VER}.tar.gz && \
   mv -v raven-${RAVEN_VER} /raven && \
   mkdir /data

# Build the Raven release and test (does NOT build python bindings)
RUN cd /raven && \
   cmake -S ./ -B./build \
    -DRAVEN_BUILD_EXE=1 \
    -DCMAKE_BUILD_TYPE=Release \
    -DRAVEN_BUILD_TESTS=1 && \
   cmake --build ./build && \
   cmake --install ./build

# Use for e.g. $PATH and locale settings for compatibility with Singularity
ENV PATH="/raven-${SOFTWARENAME_VER}/bin:$PATH" \
 LC_ALL=C

RUN raven --version

WORKDIR /data