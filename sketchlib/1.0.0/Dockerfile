FROM mambaorg/micromamba:debian13-slim

# Install sketchlib (rust binary) and python3 via conda
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        sketchlib=0.2.4 \
        python=3.12 \
    && micromamba clean --all --yes

# Ensure conda env is on PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Verify both are available
RUN sketchlib --version && python3 --version

USER $MAMBA_USER

CMD ["bash"]