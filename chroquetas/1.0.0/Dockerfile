FROM ubuntu:jammy as app

ARG CHROQUETAS_VER="1.0.0"
ARG MINIPROT_VER="0.18"
ARG SEQKIT_VER="2.10.1"

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="ChroQueTas"
LABEL software.version="${CHROQUETAS_VER}"
LABEL description="Fungal AMR detection"
LABEL website="https://github.com/theiagen_docker_builds/chroquetas"
LABEL license="https://github.com/theiagen_docker_builds/blob/main/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

# Install dependencies available in Debian repos
RUN apt-get update && apt-get install -y --no-install-recommends \
 ca-certificates \
 ncbi-blast+ \
 emboss \
 mafft \
 build-essential \
 bzip2 \
 wget && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# install ChroQueTas 
RUN wget https://github.com/nmquijada/ChroQueTas/archive/refs/tags/v${CHROQUETAS_VER}.tar.gz && \
    tar -xzf v${CHROQUETAS_VER}.tar.gz && \
    rm v${CHROQUETAS_VER}.tar.gz && \
    mv ChroQueTas-${CHROQUETAS_VER} /chroquetas && \
    chmod +x /chroquetas/bin/ChroQueTas.sh

# install FungAMR database - has separate versions, but tied to ChroQueTas release, so we glob
RUN tar -xvzf /chroquetas/FungAMR_db/db*.tgz -C /chroquetas/FungAMR_db/ && \
    rm /chroquetas/FungAMR_db/db*.tgz && \
    mv /chroquetas/FungAMR_db/db*/* /chroquetas/FungAMR_db/ && \
    rmdir /chroquetas/FungAMR_db/db*

# install ChroQueTas dependencies not available in Debian repos
# miniprot
RUN wget https://github.com/lh3/miniprot/releases/download/v0.18/miniprot-${MINIPROT_VER}_x64-linux.tar.bz2 && \
 tar -xjf miniprot-${MINIPROT_VER}_x64-linux.tar.bz2 && \
 rm miniprot-${MINIPROT_VER}_x64-linux.tar.bz2 && \
 ln -s /miniprot-${MINIPROT_VER}_x64-linux/miniprot /chroquetas/bin/miniprot

# seqkit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VER}/seqkit_linux_amd64.tar.gz && \
    tar -xvzf seqkit_linux_amd64.tar.gz && \
    rm seqkit_linux_amd64.tar.gz && \
    mv seqkit /chroquetas/bin/seqkit

ENV PATH="/chroquetas/bin:$PATH" \
 LC_ALL=C

CMD [ "ChroQueTas.sh" ]

RUN mkdir /data
WORKDIR /data

# run a test
FROM app as test

WORKDIR /test

# Candida albicans test
RUN ChroQueTas.sh -g /chroquetas/test_dataset/Calbicans_SRR13587609.fasta.gz -s Candida_albicans -o /test/Calbicans_SRR13587609_ChroQueTas && \
    rm -rf /test