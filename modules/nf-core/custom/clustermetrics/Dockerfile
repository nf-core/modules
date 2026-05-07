FROM nfcore/base:2.0

LABEL authors="dbaku42" \
      description="Docker image containing all requirements for nf-core/custom/clustermetrics"

COPY environment.yml /
RUN micromamba install -y -n base -f /environment.yml && \
    micromamba clean -a -y
