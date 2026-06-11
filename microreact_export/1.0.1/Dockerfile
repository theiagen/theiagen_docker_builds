FROM python:3.13-slim AS app

ARG PANDAS_VER="2.3.3"
ARG REQUESTS_VER="2.32.5"

LABEL description="Docker container for Microreact export script"
LABEL website="https://github.com/theiagen_docker_builds/microreact_export"
LABEL license="https://github.com/theiagen_docker_builds/blob/main/LICENSE"
LABEL maintainer="Andrew Hale"
LABEL maintainer.email="andrew.hale@theiagen.com"
LABEL dockerfile.version="1.0.1"

RUN apt-get update && apt-get install -y --no-install-recommends gcc\
    && pip3 install --no-cache-dir \
        pandas==${PANDAS_VER} \
        requests==${REQUESTS_VER} \
    && apt-get remove -y gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /scripts

COPY microreact_export.py /scripts/microreact_export.py

ENTRYPOINT [ "python", "/scripts/microreact_export.py" ]
