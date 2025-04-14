FROM ubuntu:jammy AS app

ARG MEGAHIT_VER="1.2.9"

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="MEGAHIT"
LABEL software.version="${MEGAHIT_VER}"
LABEL description="de novo metagenome assembler"
LABEL website="https://github.com/voutcn/megahit"
LABEL license="https://github.com/voutcn/megahit/main/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

# install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends -y \
    python3 \
    bzip2 \
    wget \
    ca-certificates \
    gzip && \
    ln -s /usr/bin/python3 /usr/bin/python

# download latest releast
RUN wget -q https://github.com/voutcn/megahit/releases/download/v${MEGAHIT_VER}/MEGAHIT-${MEGAHIT_VER}-Linux-x86_64-static.tar.gz && \
    tar -xzf MEGAHIT-${MEGAHIT_VER}-Linux-x86_64-static.tar.gz && \
    rm MEGAHIT-${MEGAHIT_VER}-Linux-x86_64-static.tar.gz && \
    mkdir /data

# set PATH and locale settings for singularity
ENV LC_ALL=C.UTF-8 \
  PATH="${PATH}:/MEGAHIT-${MEGAHIT_VER}-Linux-x86_64-static/bin"

CMD ["megahit", "--help"]

WORKDIR /data

RUN megahit --test && megahit --test --kmin-1pass
ENTRYPOINT ["megahit"]