FROM python:3.11-slim AS app

WORKDIR /data

# 'LABEL' instructions tag the image with metadata that might be important to the user
LABEL base.image="python:3.11-slim"
LABEL dockerfile.version="3"
LABEL software="AlleleClustering"
LABEL description="PulseNet 2.0 Allele Clustering algorithm: Distance matrix calculation and tree inference code for PulseNet 2.0 application"
LABEL website="https://github.com/ncezid-biome/pulsenet2.0-trees"

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps && \
  apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir numpy==2.4.4

COPY src/*.py /data/

# 'CMD' instructions set a default command when the container is run. This is typically 'tool --help.'
CMD [ "python3", "AlleleClustering.py", "--help" ]
