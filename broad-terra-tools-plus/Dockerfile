# Dockerfile modified from here: https://github.com/broadinstitute/terra-tools/blob/master/Dockerfile
# updated 2023-02-13 to bring more tool up to date

FROM google/cloud-sdk:417.0.1
	
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

#### not sure why python was installed manually, leaving cruft here for now  ######
    # && wget https://www.python.org/ftp/python/3.7.4/Python-3.7.4.tgz \
    # && tar xzf Python-3.7.4.tgz \
    # && cd Python-3.7.4 \
    # && ./configure --enable-optimizations \
    # && make altinstall \
    # && update-alternatives --install /usr/bin/python3 python3 /usr/bin/Python-3.7.4/python 1 

# install pip3   
# RUN rm /usr/bin/lsb_release \ 
#     && apt-get update -y \
#     && apt-get install python3-pip idle3 -y \ 
#     && python3 -m pip install --no-cache-dir --upgrade pip

#COPY scripts /scripts

# clone broad/terra-tools repo
RUN git clone https://github.com/broadinstitute/terra-tools.git && \
 mkdir /scripts && \
 cp -vr /terra-tools/scripts /scripts

# original requirements.txt:
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

#COPY requirements.txt .
#RUN pip3 install -r requirements.txt

#ENV PYTHONPATH "/scripts:${PYTHONPATH}"

#CMD ["/bin/bash"]