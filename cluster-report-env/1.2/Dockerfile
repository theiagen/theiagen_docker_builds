# base image
FROM ubuntu:xenial

# metadata
LABEL base.image="ubuntu:xenial"
LABEL description="An r environment and script for generating a genome cluster report"
LABEL website="https://github.com/kevinlibuit/reports"
LABEL comments: "Include template rmd file; Dockerfile based on https://github.com/StaPH-B/docker-builds/blob/master/cluster-report-env/1.1/Dockerfile"
LABEL maintainer="Kevin Libuit"
LABEL maintainer.email="Kevin@libuitsci.com"

# install ubuntu dependencies
RUN apt-get update && \
  apt-get -y install software-properties-common \
  apt-transport-https \
  ca-certificates \
  gnupg \
  curl \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2 \
  libxml2-dev \
  pandoc \
  pandoc-citeproc \
  wget \
  git \
  imagemagick \
  libmagick++-dev \
  texlive-latex-base \
  texlive-fonts-recommended \
  texlive-latex-recommended \
  texlive-latex-extra \
  texlive-xetex \
  fonts-roboto

# add keys and ppa; update sources; install latest version of R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
  apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
  add-apt-repository ppa:marutter/c2d4u3.5 && \
  echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list && \
  apt-get update && apt-get -y upgrade && apt-get -y install r-base \
  r-base-dev

# install R packages
RUN R -e "install.packages(c('knitr',\
  'rmarkdown',\
  'tidyverse',\
  'devtools',\
  'phytools',\
  'viridisLite',\
  'viridis',\
  'kableExtra',\
  'pander'), repos = 'http://cran.us.r-project.org')" && \
  R -e "devtools::install_github('AbigailShockey/ggtree')" && \
  R -e "tinytex::install_tinytex()"

# install report repository & grab template file report
RUN wget https://github.com/theiagen/reports/archive/v1.2.1.tar.gz && \
  tar -xvf v1.2.1.tar.gz && \
  rm -rf v1.2.1.tar.gz && \
  mkdir /reports && \
  cp reports-1.2.1/*report_template.Rmd /reports

WORKDIR /data
ENV PATH "/reports:$PATH"
