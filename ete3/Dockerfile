FROM python:3.9.17-slim as app

# List all software versions are ARGs near the top of the dockerfile
ARG ETE3="3.1.3"

# 'LABEL' instructions tag the image with metadata that might be important to the user
LABEL base.image="python:3.9.17-slim"
LABEL dockerfile.version="1"
LABEL software="ETE3"
LABEL software.version="${ETE3}"
LABEL website="https://github.com/theiagen/theiagen_docker_builds/"
LABEL license="https://github.com/theiagen/theiagen_docker_builds/blob/master/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

# Install and/or setup more things. Make /data for use as a working dir
# For readability, limit one install per 'RUN' statement.
RUN pip install --no-cache ete3==${ETE3} six numpy

# 'ENV' instructions set environment variables that persist from the build into the resulting image
ENV PATH="/software-${ETE3}/bin:$PATH" \
    LC_ALL=C

# 'CMD' instructions set a default command when the container is run. This is typically 'tool --help.'
CMD [ "ete3", "--help" ]

# 'WORKDIR' sets working directory
WORKDIR /data