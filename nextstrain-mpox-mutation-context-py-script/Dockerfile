# base image
FROM ubuntu:jammy

# metadata
LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="mutation_context.py from Nextstrain"
LABEL software.version="none"
LABEL description="a docker image for running mutation_context.py from Augur Mpox public build"
LABEL website="https://github.com/nextstrain/mpox/blob/master/nextclade/scripts/mutation_context.py"
LABEL license="https://github.com/nextstrain/mpox/blob/master/LICENSE"
LABEL maintainer="Curtis Kapsak"
LABEL maintainer.email="curtis.kapsak@theiagen.com"

# install dependencies
RUN apt-get update && apt-get -y --no-install-recommends install \
    python3 \
    wget \
    ca-certificates && \
    rm -rf /var/lib/apt/lists/* && apt-get autoclean

# install the mutation_context.py script to be located at /scripts/mutation_context.py
RUN mkdir /scripts && cd /scripts && \
  wget https://raw.githubusercontent.com/nextstrain/mpox/master/nextclade/scripts/mutation_context.py && \
  chmod +x /scripts/mutation_context.py

# set paths & locale for singularity compatibilty
ENV PATH="/scripts:\
${PATH}" \
LC_ALL=C

# set final working directory as /data
WORKDIR /data