FROM mambaorg/micromamba:1.5.8 AS app

ARG AMRRULES_BRANCH="genome_summary_report_dev"
ARG AMRRULES_VERSION="beta"

LABEL base.image="mambaorg/micromamba:1.5.8"
LABEL software="AMRrules"
LABEL software.version="${AMRRULES_VERSION}"
LABEL description="AMRrules: interpretation rules for AMRFinderPlus results"
LABEL website="https://github.com/AMRverse/AMRrules"
LABEL license="https://github.com/AMRverse/AMRrules/blob/main/LICENSE"
LABEL maintainer="Michal Babinski"
LABEL maintainer.email="michal.babinski@theiagen.com"

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    make \
    ca-certificates \
    && apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN micromamba create -y -n amrrules_beta -c bioconda -c conda-forge python=3.12 pip && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/envs/amrrules_beta/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=amrrules_beta

WORKDIR /opt
RUN git clone https://github.com/AMRverse/AMRrules && \
    cd AMRrules && \
    git checkout ${AMRRULES_BRANCH} && \
    make dev

RUN amrrules --download-resources

WORKDIR /data

FROM app AS test

RUN amrrules --help
