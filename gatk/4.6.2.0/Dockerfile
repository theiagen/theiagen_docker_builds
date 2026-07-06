FROM mambaorg/micromamba AS app

LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"
LABEL description="GATK Variant Calling Pipeline"
LABEL version="4.6.2.0"

# copy version constraints
COPY environment.yml /tmp/environment.yml

# install GATK
RUN micromamba install -y -f /tmp/environment.yml && \
    micromamba clean -y --all

# setup run environment
ENV PATH="/opt/conda/bin:$PATH" \
    CONDA_PREFIX="/opt/conda"

WORKDIR /data