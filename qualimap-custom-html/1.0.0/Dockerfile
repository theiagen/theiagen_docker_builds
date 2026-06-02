FROM continuumio/miniconda3:latest
ARG QUALIMAP_VER="2.3"


RUN conda update -n base -c defaults conda -y && \
    conda install -n base -c conda-forge mamba -y

RUN mamba install --name base -c conda-forge -c bioconda \
    bioconductor-repitools \
    qualimap=${QUALIMAP_VER} \
    pandas \
    plotly && \
    mamba clean -a -f -y

COPY plot_coverage.py /usr/local/bin/plot_coverage.py
RUN chmod +x /usr/local/bin/plot_coverage.py

ENV PATH="/opt/conda/bin/:${PATH}" \
    LC_ALL=C.UTF-8

CMD [ "/bin/bash" ]