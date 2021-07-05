FROM mambaorg/micromamba:0.14.0
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes && rm /tmp/environment.yml && \
    # Install ps to use trace in nextflow
    apt-get update && \
    apt install -y procps g++ --no-install-recommends --no-install-suggests && \
    apt-get clean && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
