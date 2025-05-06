FROM python:3.13-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /scripts

COPY assembly_shovilter.py /scripts/assembly_shovilter.py

RUN chmod +x /scripts/assembly_shovilter.py

# Only depend on biopython
RUN pip install --no-cache-dir biopython

ENTRYPOINT ["python", "/scripts/assembly_shovilter.py"]
