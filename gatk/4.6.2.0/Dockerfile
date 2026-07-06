FROM openjdk:17.0.1-jdk-slim as app

ARG GATK_VER="4.6.2.0"

# 'LABEL' instructions tag the image with metadata that might be important to the user
LABEL base.image="openjdk:17.0.1-jdk-slim"
LABEL dockerfile.version="1"
LABEL software="GATK"
LABEL software.version="${GATK_VER}"
LABEL description="GATK variant analysis and genomics"
LABEL website="https://github.com/theiagen/theiagen_docker_builds/"
LABEL license="https://github.com/theiagen/theiagen_docker_builds/blob/master/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

ARG DEBIAN_FRONTEND=noninteractive

# Install dependencies via apt-get
RUN apt update && \
  apt install -y --no-install-recommends \
  wget \
  git \
  git-lfs \
  python \
  r-base && \
  rm -rf /var/lib/apt/lists/*

# Download GATK - must be cloned to build
RUN git clone --branch ${GATK_VER} https://github.com/broadinstitute/gatk/

# Build GATK only for local run
RUN cd gatk && \
  export GRADLE_OPTS="-Xmx4048m -Dorg.gradle.daemon=false" && \
  ./gradlew localJar 

# 'CMD' instructions set a default command when the container is run. This is typically 'tool --help.'
CMD [ "gatk", "--help" ]

# 'WORKDIR' sets working directory
WORKDIR /data

FROM app as test

# set working directory so that all test inputs & outputs are kept in /test
WORKDIR /test