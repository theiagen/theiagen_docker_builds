FROM python:3.11-slim AS app

# 'LABEL' instructions tag the image with metadata that might be important to the user
LABEL base.image="python:3.11-slim"
LABEL dockerfile.version="1"
LABEL software="AlleleClustering"
LABEL description="PulseNet 2.0 Allele Clustering algorithm: Distance matrix calculation and tree inference code for PulseNet 2.0 application"
LABEL website="https://github.com/ncezid-biome/pulsenet2.0-trees"

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps && \
  apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache numpy>=2.4.4

COPY *.py /data

# 'ENV' instructions set environment variables that persist from the build into the resulting image
# Use for e.g. $PATH and locale settings for compatibility with Singularity
ENV PATH="/software-${SOFTWARENAME_VER}/bin:$PATH" \
    LC_ALL=C

# 'CMD' instructions set a default command when the container is run. This is typically 'tool --help.'
CMD [ "tool", "--help" ]

# 'WORKDIR' sets working directory
WORKDIR /data

##### ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #####
##### Step 2. Set up the testing stage.                                 #####
##### The docker image is built to the 'test' stage before merging, but #####
##### the test stage (or any stage after 'app') will be lost.           #####
##### ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #####

# A second FROM insruction creates a new stage
FROM app as test

# set working directory so that all test inputs & outputs are kept in /test
WORKDIR /test

# print help and version info; check dependencies (not all software has these options available)
# Mostly this ensures the tool of choice is in path and is executable
RUN softwarename --help && \
    softwarename --check && \
    softwarename --version

# Demonstrate that the program is successfully installed - which is highly dependant on what the tool is.

# Run the program's internal tests if available, for example with SPAdes:
RUN spades.py --test

# Option 1: write your own tests in a bash script in the same directory as your Dockerfile and copy them:
COPY my_tests.sh .
RUN bash my_tests.sh

# Option 2: write below common usage cases, for example with tb-profiler:
RUN wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_1.fastq.gz && \
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_2.fastq.gz && \
    tb-profiler profile -1 ERR1664619_1.fastq.gz -2 ERR1664619_2.fastq.gz -t 4 -p ERR1664619 --txt
