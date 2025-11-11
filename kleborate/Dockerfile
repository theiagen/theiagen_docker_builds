FROM continuumio/miniconda3:latest

# Set metadata
LABEL maintainer="michal.babinski@theiagen.com"
LABEL description="Docker image for Kleborate - a tool for Klebsiella genome characterization"
LABEL version="3.2.4"

# Set working directory
WORKDIR /data

# Update conda and install mamba for faster dependency resolution
RUN conda update -n base -c defaults conda -y && \
    conda install -n base -c conda-forge mamba -y

# Create conda environment with all dependencies
RUN mamba create -n klebsiella_analysis \
    -c bioconda -c conda-forge \
    python=3.9 \
    minimap2 \
    mash \
    ezclermont \
    ectyper \
    stxtyper \
    ncbi-amrfinderplus \
    kleborate \
    -y

# Initialize conda for bash shell
RUN conda init bash

# Download AMRFinder database
RUN /bin/bash -c "source activate klebsiella_analysis && \
    amrfinder --update"

# Set up environment activation in bashrc
RUN echo "source activate klebsiella_analysis" >> ~/.bashrc

# Make conda environment activate on container start
ENV PATH=/opt/conda/envs/klebsiella_analysis/bin:$PATH

# Copy the parsing script into the container
COPY parse_kleborate_kleb.py /usr/local/bin/parse_kleborate_kleb.py
RUN chmod +x /usr/local/bin/parse_kleborate_kleb.py

# Set the default command to bash
CMD ["/bin/bash"]

# Verify installation (optional - can be removed for smaller image)
RUN /bin/bash -c "source activate klebsiella_analysis && \
    kleborate --version && \
    amrfinder --version"