FROM ubuntu:jammy AS app

ARG ARLN_STATS_VERSION="0.1.0"

# Metadata
LABEL base.image="ubuntu:jammy"
LABEL software="ARLN_Stats"
LABEL software.version="${ARLN_STATS_VERSION}"
LABEL description="Docker image for task_arln_stats.wdl used to agregate and calculate statistics needed for ARLN submissions"
LABEL maintainer="Andrew Hale"
LABEL maintainer.email="andrew.hale@theiagen.com"

# Install Python and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-venv \
    python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/bin/python3 /usr/bin/python \
    && ln -sf /usr/bin/pip3 /usr/bin/pip

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libz-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for FASTQ handling
RUN pip install --no-cache-dir \
    biopython

# Create and set working directory
WORKDIR /data
COPY *.fastq.gz NCBI_Assembly_stats_20240124.txt ./

WORKDIR /scripts
COPY *.py ./

# Step 2: Testing stage
FROM app AS test

WORKDIR /test

RUN mv ../data/*.fastq.gz /test/ && \
    python ../scripts/q30.py -i *_R1_001.fastq.gz >> q30_test.txt && \
    python ../scripts/q30.py -i *_R1_001.fastq.gz >> q30_test.txt && \
    cat q30_test.txt