FROM staphb/clair3:1.0.10

USER root
WORKDIR /

RUN apt-get update && apt-get install -y \
    git \
    python3 \
    && rm -rf /var/lib/apt/lists/*

# Get latest models from rerio repo
WORKDIR /rerio
RUN git clone https://github.com/nanoporetech/rerio . && \
    for model in \
        "clair3_models/r1041_e82_400bps_sup_v500_model" \
        "clair3_models/r1041_e82_400bps_hac_v500_model" \
        "clair3_models/r1041_e82_400bps_sup_v410_model" \
        "clair3_models/r1041_e82_400bps_hac_v410_model"; \
    do \
        python3 download_model.py --clair3 $model; \
    done && \
    cp -r clair3_models/* /clair3/models/ && \
    ls -la /clair3/models/

WORKDIR /data