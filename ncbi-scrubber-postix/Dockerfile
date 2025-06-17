FROM ncbi/sra-human-scrubber:latest AS build

# Install additional system utilities for debugging/monitoring
RUN apt-get update && apt-get install -y \
    procps \
    htop \
    tree \
    curl \
    wget \
    vim \
    nano \
    google-cloud-sdk \
    && rm -rf /var/lib/apt/lists/*

COPY . /opt/scrubber/
WORKDIR /opt/scrubber
RUN /opt/scrubber/init_db.sh
ENV PATH="/usr/bin:/opt/scrubber/bin:/opt/scrubber/scripts:${PATH}"
RUN ldconfig
RUN ./scripts/scrub.sh -p 6 -t