# building from staphb/irma as the base
FROM staphb/irma:1.3.1 

# Metadata
LABEL base.image="staphb/irma:1.3.1"
LABEL software="IRMA"
LABEL software.version="1.3.1"
LABEL maintainer="Andrew Hale"
LABEL maintainer.email="andrew.hale@theiagen.com"

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip

WORKDIR /scripts

COPY irma_helper.py /scripts/irma_helper.py

RUN pip install --no-cache-dir biopython pandas

RUN chmod +x /scripts/irma_helper.py