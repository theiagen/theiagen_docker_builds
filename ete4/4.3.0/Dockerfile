# 'FROM' defines the base docker image. This command has to come first in the file
# The 'as' keyword lets you name the folowing stage. The production image uses everything to the 'app' stage.
FROM mambaorg/micromamba:2.3.0-ubuntu22.04 AS app

# List all software versions are ARGs near the top of the dockerfile
# 'ARG' sets environment variables during the build stage
# 'ARG' variables are ONLY available during image build, they do not persist in the final image
ARG ETE4_VERSION="4.3.0"

# build and run as root users since micromamba image has 'mambauser' set as the $USER
USER root
# set workdir to default for building; set to /data at the end
WORKDIR /

# 'LABEL' instructions tag the image with metadata that might be important to the user

LABEL base.image="mambaorg/micromamba:2.3.0-ubuntu22.04"
LABEL dockerfile.version="1"
LABEL software="ETE4"
LABEL software.version="${ETE4_VERSION}"
LABEL description="Hosts the ETE toolkit"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

# Install your desired software into the base conda/micromamba environment, pinning the version
RUN micromamba install --name base conda-forge::ete4=${ETE4_VERSION} && \
 micromamba clean -a -f -y && \
 mkdir /data

# 'ENV' instructions set environment variables that persist from the build into the resulting image
# set the environment, add base conda/micromamba bin directory into path
# set locale settings to UTF-8
ENV PATH="/opt/conda/bin/:${PATH}" \
 LC_ALL=C.UTF-8

# Pre-download the NCBI taxonomy database to avoid first-run download delays
RUN python -c "from ete4 import NCBITaxa; ncbi = NCBITaxa()"

CMD [ "ete4", "--help" ]

# set final working directory to /data
WORKDIR /data