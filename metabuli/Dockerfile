FROM ubuntu:jammy as app

ARG METABULI_VER="1.1.1"

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="Metabuli"
LABEL software.version="${METABULI_VER}"
LABEL description="Taxonomic read classification"
LABEL website="https://github.com/theiagen_docker_builds/metabuli"
LABEL license="https://github.com/theiagen_docker_builds/blob/main/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

# Install dependencies available in Debian repos
RUN apt-get update && apt-get install -y --no-install-recommends \
 ca-certificates \
 build-essential \
 curl \
 xxd \
 cmake \
 zlib1g-dev \
 libbz2-dev \
 libatomic1 \
 git \
 wget && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Download Metabuli
RUN git clone --recurse-submodules -b ${METABULI_VER} https://github.com/steineggerlab/Metabuli.git /metabuli


# Build Metabuli
RUN cd /metabuli; \
  mkdir build && \
  cd build && \
  cmake -DCMAKE_BUILD_TYPE=Release -DNATIVE_ARCH=0 -DHAVE_AVX2_1=1 .. && \
  make -j $(nproc)

ENV PATH="/metabuli/build/src:$PATH" \
 LC_ALL=C

CMD [ "metabuli -h" ]

RUN mkdir /data
WORKDIR /data