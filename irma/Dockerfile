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
    python3