FROM ubuntu:jammy AS app

# VEP_VER is the git release tag (release/${VEP_VER}); it is also the software.version
ARG VEP_VER="116.0"
# SAMTOOLS_VER is the samtools release tag; compiled from source below
ARG SAMTOOLS_VER="1.23"
ARG THEIAGENE_VER="1.0.1"

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="ensembl-vep"
LABEL software.version="${VEP_VER}"
LABEL description="Ensembl Variant Effect Predictor (VEP): determines the effect of variants on genes, transcripts, protein sequence, and regulatory regions"
LABEL website="https://github.com/Ensembl/ensembl-vep"
LABEL license="https://github.com/Ensembl/ensembl-vep/blob/main/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

ARG DEBIAN_FRONTEND=noninteractive

# install build + runtime dependencies and cleanup apt garbage
# - build-essential/autoconf/automake/libperl-dev: compile htslib + Bio::DB::HTS
# - zlib/bz2/lzma/curl/ssl -dev: required to build htslib
# - Archive::Zip + DBI: required by INSTALL.pl before running
# - DBD::mysql: database/cache lookups; Set::IntervalTree/JSON/PerlIO::gzip: VEP performance + output features
# - Module::Build: builds VEP's XS Perl module dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    ca-certificates \
    git \
    wget \
    curl \
    unzip \
    python3 \
    python3-pip \
    perl \
    libperl-dev \
    cpanminus \
    build-essential \
    autoconf \
    automake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses-dev \
    libmodule-build-perl \
    libarchive-zip-perl \
    libdbi-perl \
    libdbd-mysql-perl \
    libset-intervaltree-perl \
    libjson-perl \
    libperlio-gzip-perl \
    liblist-moreutils-perl \
    libwww-perl \
    libtry-tiny-perl && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# jammy ships pip 22.0.2 / setuptools 59.6.0, which predate PEP 621 support and
# silently build an empty "UNKNOWN-0.0.0" wheel from theiagene's [project] table
# (import theiagene then fails). Upgrade pip/setuptools/wheel first so the
# pyproject.toml metadata is read correctly.
RUN python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel

# clone ensembl-vep, pin to the release/${VEP_VER} tag, and install the API + Faidx/htslib
# --AUTO a : install Ensembl API and Faidx/htslib   --NO_UPDATE : do not check for updates
# --NO_TEST : skip the API test suite (requires network access to Ensembl during build)
# caches, FASTA files, and plugins are intentionally NOT installed here; mount/download them at runtime
RUN mkdir -p /opt/vep/src && cd /opt/vep/src && \
    git clone https://github.com/Ensembl/ensembl-vep.git && \
    cd ensembl-vep && \
    git checkout release/${VEP_VER} && \
    perl INSTALL.pl --AUTO a --NO_UPDATE --NO_TEST && \
    rm -rf /opt/vep/src/ensembl-vep/.git

# download, compile, and install samtools (bundles its own htslib in the release tarball)
# samtools is handy for prepping VEP inputs (indexing/subsetting BAM/CRAM, faidx on reference FASTA)
RUN cd /tmp && \
    wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER} && \
    ./configure && \
    make && \
    make install && \
    cd /tmp && rm -rf samtools-${SAMTOOLS_VER} samtools-${SAMTOOLS_VER}.tar.bz2

# download theiagene release and install it into the system python via pip
RUN cd /tmp && \
    wget -q https://github.com/theiagen/theiagene/archive/refs/tags/v${THEIAGENE_VER}.tar.gz && \
    tar -xzf v${THEIAGENE_VER}.tar.gz && \
    cd theiagene-${THEIAGENE_VER} && \
    python3 -m pip install --no-cache-dir . && \
    cd /tmp && rm -rf theiagene-${THEIAGENE_VER} v${THEIAGENE_VER}.tar.gz

# VEP executables (vep, filter_vep, haplo, variant_recoder) and the htslib binaries
# (tabix, bgzip) that INSTALL.pl compiled into ./htslib
ENV PATH="/opt/vep/src/ensembl-vep:/opt/vep/src/ensembl-vep/htslib:${PATH}" \
    PERL5LIB="/opt/vep/src/ensembl-vep:${PERL5LIB}" \
    LC_ALL=C

# final working directory is /data
WORKDIR /data

# default command is to pull up help options
CMD ["vep", "--help"]

### start of test stage ###
FROM app AS test

# confirm the VEP executables installed and are on PATH
RUN vep --help && \
    filter_vep --help && \
    haplo --help && \
    variant_recoder --help

# confirm htslib tooling compiled and is on PATH
RUN tabix --version && bgzip --version

# confirm samtools compiled and is on PATH
RUN samtools --version

# confirm theiagene imports and its CLI entrypoint is on PATH
RUN python3 -c "import theiagene" && \
    theiagene --help

# confirm the Ensembl API and key Perl modules load (Bio::DB::HTS is compiled by INSTALL.pl)
# and report the installed API version deterministically
RUN perl -e 'use Bio::EnsEMBL::Registry; use Bio::EnsEMBL::Variation::Variation; use Bio::EnsEMBL::IO::Parser; use Bio::DB::HTS; use DBI; use DBD::mysql; use Archive::Zip; use Set::IntervalTree; use JSON; use PerlIO::gzip; print "all VEP perl dependencies load OK\n"; print "Ensembl API version: ", Bio::EnsEMBL::Registry->software_version, "\n"'
