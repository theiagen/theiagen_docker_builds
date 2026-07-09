FROM mambaorg/micromamba AS app

LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"
LABEL description="GATK Variant Calling Pipeline"
LABEL version="4.6.2.0"

# copy version constraints
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# install GATK
RUN micromamba install -y -f /tmp/environment.yml && \
    micromamba clean -y --all && \
    rm /tmp/environment.yml

# setup run environment
ENV PATH="/opt/conda/bin:$PATH" \
    CONDA_PREFIX="/opt/conda"

USER root

RUN mkdir /data && \
  chown $MAMBA_USER:$MAMBA_USER /data

USER $MAMBA_USER

WORKDIR /data
