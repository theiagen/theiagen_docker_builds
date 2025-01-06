FROM mambaorg/micromamba:2.0.5

LABEL maintainer="Michal Babinski"
LABEL maintainer.email="michal.babinski@theiagen.com"
LABEL description="ARTIC Pipeline for nanopore viral genome sequencing"
LABEL version="1.5.8"

# Install ARTIC pipeline
RUN micromamba install -y -c bioconda -c conda-forge \
    artic \
    git \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:$PATH"

WORKDIR /data

RUN git clone --recursive https://github.com/artic-network/primer-schemes