FROM mambaorg/micromamba AS app

LABEL maintainer="Michal Babinski"
LABEL maintainer.email="michal.babinski@theiagen.com"
LABEL description="ARTIC Pipeline for nanopore viral genome sequencing"
LABEL version="1.9.0"

USER root

# Install ARTIC pipeline
COPY environment.yml /tmp/environment.yml

RUN micromamba install -y -f /tmp/environment.yml && \
    micromamba clean -y --all && \
    rm -f /tmp/environment.yml

ENV PATH="/opt/conda/bin:$PATH" \
    CONDA_PREFIX="/opt/conda"

# Get models for newer flowcell chemistries
RUN artic_get_models

WORKDIR /data