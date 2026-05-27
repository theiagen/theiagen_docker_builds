FROM staphb/mlst:2.23.0-2024-12-31

ARG DATABASE_VERSION="2.33.1"
ARG MLST_VERSION="2.23.0"

# Metadata
LABEL base.image="staphb/mlst:2.23.0-2024-12-31"
LABEL software="MLST"
LABEL software.version="${MLST_VERSION}"
LABEL database.version="${DATABASE_VERSION}}"
LABEL maintainer="Andrew Hale"
LABEL maintainer.email="andrew.hale@theiagen.com"

RUN apt-get update && apt-get -y --no-install-recommends install \
    unzip

WORKDIR /work

# Pull 2.33.1 repo to get downgraded db
RUN wget https://github.com/tseemann/mlst/archive/refs/tags/v${DATABASE_VERSION}.zip
RUN unzip v${DATABASE_VERSION}.zip

# Remove non compliant db and replace with downgraded db
RUN rm -rf ../mlst-2.23.0/db/
RUN cp -r mlst-2.33.1/db ../mlst-2.23.0/ 

# Remove downloaded repo
RUN rm -rf mlst-2.33.1

WORKDIR /data

# Test on bundled FASTA
RUN mlst ../mlst-2.23.0/test/example.fna.gz 
# This file is present only in the downgraded db not in 2.23.0; this verifies we have the correct db
RUN cat ../mlst-2.23.0/db/VERSION.txt