FROM mambaorg/micromamba:0.14.0
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes && rm /tmp/environment.yml
