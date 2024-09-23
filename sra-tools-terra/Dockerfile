FROM google/cloud-sdk:427.0.0
	
ARG VERSION="current"
# "current" is the latest version, not the actual version number. this is very smart

LABEL software="sra-tools"
LABEL software.version="${VERSION}"
LABEL description="Downloads biological sequence data from NCBI"
LABEL website="https://github.com/ncbi/datasets"
LABEL documentation="https://www.ncbi.nlm.nih.gov/datasets/docs/v1/"
LABEL license="https://github.com/ncbi/datasets/blob/master/pkgs/ncbi-datasets-cli/LICENSE.md"

# unzip isn't needed for datasets/dataformat, but it is often used after downloading files with datasets
RUN apt-get update && apt-get install -y --no-install-recommends \
  wget \
  ca-certificates \
  unzip \
  curl \
  uuid-runtime \
  libxml2-utils \
  wget \
  python3 \
  python3-pip \
  jq \
  python3-setuptools && \
  apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-centos_linux64-cloud.tar.gz | tar xz -C /
ENV PATH=/usr/local/ncbi/sra-tools/bin:${PATH}
RUN mkdir -p /root/.ncbi && \
    printf '/LIBS/IMAGE_GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg && \
    printf '/libs/cloud/report_instance_identity = "true"\n' >> /root/.ncbi/user-settings.mkfg

RUN curl https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz --output edirect.tar.gz && \
    tar xzf edirect.tar.gz && \
    rm edirect.tar.gz 


ENV PATH=/edirect:${PATH}

# Tell gcloud to save state in /.config so it's easy to override as a mounted volume.
ENV HOME=/

WORKDIR /

# clone broad/terra-tools repo
# copying scripts dir to `/scripts` for consistency
RUN git clone https://github.com/broadinstitute/terra-tools.git && \
  mkdir /scripts && \
  cp -vr /terra-tools/scripts/* /scripts

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
