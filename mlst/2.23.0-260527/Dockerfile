FROM staphb/mlst:2.23.0-2024-12-31

# Metadata
LABEL base.image="staphb/mlst:2.23.0-2024-12-31"
LABEL software="MLST"
LABEL software.version="2.23.0"
LABEL maintainer="Andrew Hale"
LABEL maintainer.email="andrew.hale@theiagen.com"

RUN rm -rf /mlst-2.23.0/db/

COPY db/ /mlst-2.23.0/db

WORKDIR /data

# Test on bundled FASTA
RUN mlst ../mlst-2.23.0/test/example.fna.gz