FROM ubuntu:jammy AS app

# ARGs for software versions
ARG MAVEN_VERSION=3.9.4
ARG AMRSEARCH_VERSION="latest"
ARG USER_HOME_DIR="/root"
ARG SHA=deaa39e16b2cf20f8cd7d232a1306344f04020e1f0fb28d35492606f647a60fe729cc40d3cba33e093a17aed41bd161fe1240556d0f1b80e773abd408686217e
ARG BASE_URL=https://apache.osuosl.org/maven/maven-3/${MAVEN_VERSION}/binaries

# Metadata
LABEL base.image="ubuntu:jammy"
LABEL software="AMRsearch"
LABEL software.version="${AMRSEARCH_VERSION}"
LABEL description="AMRsearch pipeline for antimicrobial resistance prediction"
LABEL website="https://github.com/pathogenwatch-oss/amr-search"
LABEL license="https://github.com/pathogenwatch-oss/amr-search/blob/main/LICENSE"
LABEL maintainer="Andrew Hale"
LABEL maintainer.email="andrew.hale@theiagen.com"

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    unzip \
    openjdk-17-jdk-headless \
    curl \
    python3 \
    python3-pip && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Install Python libraries
RUN pip3 install --no-cache-dir \
    matplotlib \
    pandas

# Install Maven
RUN mkdir -p /usr/share/maven /usr/share/maven/ref && \
    curl -fsSL -o /tmp/apache-maven.tar.gz ${BASE_URL}/apache-maven-${MAVEN_VERSION}-bin.tar.gz && \
    echo "${SHA}  /tmp/apache-maven.tar.gz" | sha512sum -c - && \
    tar -xzf /tmp/apache-maven.tar.gz -C /usr/share/maven --strip-components=1 && \
    rm -f /tmp/apache-maven.tar.gz && \
    ln -s /usr/share/maven/bin/mvn /usr/bin/mvn

# Set environment variables
ENV MAVEN_HOME=/usr/share/maven
ENV MAVEN_CONFIG="${USER_HOME_DIR}/.m2"

# Install BLAST
RUN mkdir -p /opt/blast && \
    curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz | \
    tar -zxC /opt/blast --strip-components=1

ENV PATH="/opt/blast/bin:$PATH"


# Obtain and build AMRsearch
RUN wget https://github.com/pathogenwatch-oss/amr-search/archive/refs/heads/main.zip && \
    unzip main.zip && \
    rm main.zip && \
    ls amr-search-main
# Obtain amr-libraries and place into libraries directory
RUN wget https://github.com/pathogenwatch/amr-libraries/archive/refs/tags/0.0.20.zip && \
    unzip 0.0.20.zip && \
    rm 0.0.20.zip && \
    cp amr-libraries-0.0.20/* amr-search-main/libraries/amr-libraries && \
    mv amr-libraries-0.0.20/* amr-search-main/libraries/amr-test-libraries && \
    rm -rf amr-libraries-0.0.20

RUN cd amr-search-main && \
    mkdir -p /build && mvn verify && \
    mkdir /paarsnp && \
    mv ./build/paarsnp.jar /paarsnp/paarsnp.jar && \
    mv ./build/databases /paarsnp && \
    rm -f /paarsnp/databases/*.fna && \
    mv ./resources/taxid.map /paarsnp/databases/

# Copy parse_amr_json.py
WORKDIR /scripts
COPY parse_amr_json.py . 

# Working directory and entrypoint
WORKDIR /data
ENTRYPOINT ["java", "-jar", "/paarsnp/paarsnp.jar"]

# Step 2: Testing stage
FROM app AS test

WORKDIR /test

# Test functionality with example FASTA file
RUN mv ../amr-search-main/test.fasta /test/ && \
    java -jar /paarsnp/paarsnp.jar -i /test/test.fasta -s 485
RUN python3 /scripts/parse_amr_json.py *.jsn test && \
    cat *_amr_results.csv