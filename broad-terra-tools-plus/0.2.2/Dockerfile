# Dockerfile modified from here: https://github.com/broadinstitute/terra-tools/blob/master/Dockerfile
# updated 2023-02-13 to bring more tool up to date

FROM google/cloud-sdk:422.0.0

# for easy upgrade later
ARG THEIAGEN_UTILITIES_VER="0.2.2"

# Tell gcloud to save state in /.config so it's easy to override as a mounted volume.
ENV HOME=/

# install stuff via apt; cleanup apt garbage
RUN apt-get update && apt-get install -y --no-install-recommends \
 ca-certificates \
 build-essential \
 checkinstall \
 libreadline-dev \
 libncursesw5-dev \
 libssl-dev \
 libsqlite3-dev \
 tk-dev \
 libgdbm-dev \
 libc6-dev \
 libbz2-dev \
 libffi-dev \
 zlib1g-dev \
 liblzma-dev \
 wget \
 python3 \
 python3-pip \
 python3-setuptools && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# grab theiagen/utilities code
# copy over scripts into /scripts within the container filesystem
RUN wget https://github.com/theiagen/utilities/archive/refs/tags/v${THEIAGEN_UTILITIES_VER}.tar.gz && \
 tar -xzf v${THEIAGEN_UTILITIES_VER}.tar.gz && \
 rm -v v${THEIAGEN_UTILITIES_VER}.tar.gz && \
 mkdir /scripts && \
 cp -v utilities-${THEIAGEN_UTILITIES_VER}/scripts/* /scripts

# clone broad/terra-tools repo
# copying scripts dir to `/scripts` for consistency
RUN git clone https://github.com/broadinstitute/terra-tools.git && \
 cp -vr /terra-tools/scripts/* /scripts

# original requirements.txt from broadinstitute/terra-tools repo:
# firecloud==0.16.26
# google-auth-httplib2==0.0.3
# gcs-oauth2-boto-plugin==2.5
# google-api-python-client==1.8.0
# google-cloud-storage==1.38.0
# google-cloud-bigquery==2.7.0
# pandas==1.3.0
# tqdm==4.35.0

# we are installing updated versions of each of the following via pip3/pypi:
RUN pip3 install firecloud \
google-auth-httplib2 \
gcs-oauth2-boto-plugin \
google-api-python-client \
google-cloud-storage \
google-cloud-bigquery \
pandas \
tqdm \
pyfaidx

# final working directory is /data
WORKDIR /data

# put broadinstitute/terra-tools scripts onto the PATH
ENV PATH=${PATH}:/scripts