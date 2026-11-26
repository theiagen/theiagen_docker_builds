FROM us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0 AS app

ARG FASTQ_DL_VER="2.0.4"
ARG SITE_PACKAGES="/usr/local/lib/python3.11/site-packages"
ARG ORIGINAL_MD5="0a6b79d9a477bbcaafb8a2b4496b9181"
ARG PATCHED_MD5="053c8c517dcff38a7ce08b3c46c91697"
ARG PATCH_FILE="/usr/share/doc/fastq-dl-enafree/fastq-dl-2.0.4-enafree.patch"

LABEL base.image="us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0"
LABEL dockerfile.version="1"
LABEL software="fastq-dl"
LABEL software.version="${FASTQ_DL_VER}"
LABEL description="fastq-dl 2.0.4 patched to skip ENA metadata queries when --provider sra is used"
LABEL website="https://github.com/rpetit3/fastq-dl"
LABEL license="https://github.com/rpetit3/fastq-dl/blob/master/LICENSE"

COPY fastq-dl-2.0.4-enafree.patch ${PATCH_FILE}

# fail the build if the base file or the patched file differs from the expected checksum
RUN echo "${ORIGINAL_MD5}  ${SITE_PACKAGES}/fastq_dl/fastq_dl.py" | md5sum -c - && \
    cd ${SITE_PACKAGES} && patch -p1 -i ${PATCH_FILE} && \
    echo "${PATCHED_MD5}  ${SITE_PACKAGES}/fastq_dl/fastq_dl.py" | md5sum -c - && \
    rm -f ${SITE_PACKAGES}/fastq_dl/__pycache__/fastq_dl.cpython-311.pyc && \
    python -m py_compile ${SITE_PACKAGES}/fastq_dl/fastq_dl.py

### start of test stage ###
FROM app AS test

COPY test_enafree.py /tmp/test_enafree.py

# print version and run the checks that need no network access
RUN fastq-dl --version && \
    fastq-dl --help > /dev/null && \
    python /tmp/test_enafree.py
