FROM python:3.13-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python libraries
RUN pip3 install --no-cache-dir \
    pandas \
    requests 

WORKDIR /scripts

COPY microreact_export.py /scripts/microreact_export.py

ENTRYPOINT [ "python", "/scripts/microreact_export.py" ]