FROM staphb/bbtools:39.38

ARG BBTOOLS_VER="39.38"

LABEL base.image="staphb/bbtools:39.38"
LABEL dockerfile.version="1"
LABEL software="BBTools"
LABEL software.version="${BBTOOLS_VER}"
LABEL description="A set of tools labeled as 'Bestus Bioinformaticus'"
LABEL website="https://github.com/bbushnell/BBTools"
LABEL documentation="https://bbmap.org/"
LABEL license="https://github.com/bbushnell/BBTools/blob/master/license.txt"
LABEL maintainer="Theron James"
LABEL maintainer.email="theron.james@theiagen.com"

RUN apt-get update && apt-get install -y python3 python3-pip \
    && pip3 install --no-cache-dir biopython \
    && rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/python3 /usr/bin/python

WORKDIR /data