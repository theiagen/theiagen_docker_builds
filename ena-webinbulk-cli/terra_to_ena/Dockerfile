FROM us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16

# https://linuxconfig.org/ubuntu-20-04-gpg-error-the-following-signatures-couldn-t-be-verified
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys C0BA5CE6DC6315A3 && \
    apt-get update

RUN apt-get update && apt-get install -y \
    software-properties-common \
    git \
    && rm -rf /var/lib/apt/lists/*

# Add OpenJDK repository for Java 17
RUN apt-get update && \
    apt-get install -y apt-transport-https ca-certificates wget gnupg && \
    wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add - && \
    echo "deb https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list && \
    apt-get update && \
    apt-get install -y temurin-17-jre && \
    rm -rf /var/lib/apt/lists/*

# Verify Java installation
RUN java -version

# Install Python dependencies
RUN pip3 install pandas numpy joblib

# This is where bulk_webincli expects the Webin-CLI JAR file in /
RUN echo "Downloading latest Webin-CLI..." && \
    curl -s https://api.github.com/repos/enasequence/webin-cli/releases/latest | \
    grep "browser_download_url" | cut -d : -f 2,3 | tr -d \" | wget -qi - && \
    mv webin-cli-* /scripts/webin-cli.jar && \
    chmod 554 /scripts/webin-cli.jar && \
    echo "Downloading latest Webin-CLI... [COMPLETE]"

RUN git clone https://github.com/enasequence/ena-bulk-webincli.git && \
    cp ena-bulk-webincli/bulk_webincli.py /usr/local/bin/bulk_webincli.py && \
    chmod 554 /usr/local/bin/bulk_webincli.py

COPY prepare_ena_data.py /usr/local/bin/
COPY localize_files.py /usr/local/bin/
COPY register_samples_ena.py /usr/local/bin/

RUN chmod +x /usr/local/bin/register_samples_ena.py
RUN chmod +x /usr/local/bin/prepare_ena_data.py
RUN chmod +x /usr/local/bin/localize_files.py

WORKDIR /data

ENTRYPOINT ["/bin/bash", "-c"]