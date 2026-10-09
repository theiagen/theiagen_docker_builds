# derived from the StaPH-B snpeff 5.4c Dockerfile:
# https://github.com/StaPH-B/docker-builds/tree/master/build-files/snpeff/5.4c
FROM ubuntu:noble AS app

# SNPEFF_VER is the snpEff release; the download URL uses '_' in place of '.' (5.4c -> 5_4c)
ARG SNPEFF_VER="5.4c"
ARG THEIAGENE_VER="1.0.1"

LABEL base.image="ubuntu:noble"
LABEL dockerfile.version="1"
LABEL software="SnpEff & SnpSift"
LABEL software.version="${SNPEFF_VER}"
LABEL description="Genetic variant annotation and effect prediction toolbox."
LABEL description.SnpSift="Used after SnpEff annotation to filter and manipulate annotated files."
LABEL website="https://github.com/pcingola/SnpEff"
LABEL documentation="https://pcingola.github.io/SnpEff"
LABEL license="https://github.com/pcingola/SnpEff/blob/master/LICENSE.md"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

ARG DEBIAN_FRONTEND=noninteractive

# install runtime dependencies and cleanup apt garbage
# - openjdk-21-jre-headless: runs snpEff.jar / SnpSift.jar
# - python3/python3-pip: theiagene
RUN apt-get update && apt-get install --no-install-recommends -y \
    openjdk-21-jre-headless \
    ca-certificates \
    wget \
    curl \
    unzip \
    python3 \
    python3-pip && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# download and unpack snpEff + SnpSift into /snpEff
RUN cd / && \
    SNPEFF_DOWNLOAD_VER=$(echo ${SNPEFF_VER} | tr '.' '_') && \
    wget -q https://snpeff-public.s3.amazonaws.com/versions/snpEff_v${SNPEFF_DOWNLOAD_VER}_core.zip && \
    unzip -q snpEff_v${SNPEFF_DOWNLOAD_VER}_core.zip && \
    rm snpEff_v${SNPEFF_DOWNLOAD_VER}_core.zip

# create `snpeff` and `snpsift` wrappers around the jars, and point the bundled
# helper scripts (e.g. buildDbNcbi.sh) at the absolute jar location
RUN printf '#!/bin/bash\nexec java -jar /snpEff/snpEff.jar "$@"\n' > /snpEff/snpeff && \
    printf '#!/bin/bash\nexec java -jar /snpEff/SnpSift.jar "$@"\n' > /snpEff/snpsift && \
    chmod +x /snpEff/snpeff /snpEff/snpsift /snpEff/snpEff.jar /snpEff/SnpSift.jar && \
    sed -i 's#java -jar snpEff.jar#java -jar /snpEff/snpEff.jar#g' /snpEff/scripts/*.sh && \
    chmod +x /snpEff/scripts/*.sh

# download theiagene release and install it into the system python via pip
# --break-system-packages: noble's python is marked EXTERNALLY-MANAGED (PEP 668)
RUN cd /tmp && \
    wget -q https://github.com/theiagen/theiagene/archive/refs/tags/v${THEIAGENE_VER}.tar.gz && \
    tar -xzf v${THEIAGENE_VER}.tar.gz && \
    cd theiagene-${THEIAGENE_VER} && \
    python3 -m pip install --no-cache-dir --break-system-packages . && \
    cd /tmp && rm -rf theiagene-${THEIAGENE_VER} v${THEIAGENE_VER}.tar.gz

ENV PATH="${PATH}:/snpEff:/snpEff/scripts" \
    LC_ALL=C

# final working directory is /data
WORKDIR /data

# default command is to pull up help options
CMD ["snpeff", "-help"]

### start of test stage ###
FROM app AS test

# confirm the snpEff/SnpSift wrappers are on PATH (both exit non-zero on help)
RUN snpeff -help 2>&1 || true && \
    snpsift -h 2>&1 || true && \
    snpeff -version

# confirm theiagene imports and its CLI entrypoint is on PATH
RUN python3 -c "import theiagene" && \
    theiagene --help

# download a prebuilt database and annotate the bundled example VCF
RUN snpeff download -v hg19 && \
    snpeff hg19 /snpEff/examples/test.chr22.vcf > annotated.vcf && \
    grep -v "^#" annotated.vcf | head -5 | tr "\t" "\ " | sed -e "s/.\{75\}/&\n/g"

# build a database from an NCBI GenBank accession
# from issue https://github.com/StaPH-B/docker-builds/issues/760
RUN buildDbNcbi.sh CP014866.1
