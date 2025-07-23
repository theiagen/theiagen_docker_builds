ARG FASTBAPS_VER="1.0.8"

FROM mambaorg/micromamba:2.3.0

ARG FASTBAPS_VER

# install fastbaps
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba clean -a -y 

# export the environment and fastbaps runpoint
ENV PATH="/opt/conda/bin/:${PATH}"
ENV PATH="/opt/conda/lib/R/library/fastbaps/:${PATH}"

# test the installation
RUN run_fastbaps -h
ENV FASTBAPS_VER=${FASTBAPS_VER}

WORKDIR /data