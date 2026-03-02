FROM quay.io/biocontainers/antismash:8.0.4--pyhdfd78af_0@sha256:e9f86a95ff31e975a82b4437f4d6b6b0c6d6ffb7c72bf5e41840b75d80ec263a

## Speeds up slightly builds: https://github.com/antismash/antismash/issues/837#issuecomment-2893858536
## From: https://github.com/antismash/containers/blob/main/standalone-lite/Dockerfile#L25
RUN mkdir /matplotlib && MPLCONFIGDIR=/matplotlib python3 -c "import matplotlib.pyplot as plt" && chmod -R a+rw /matplotlib

# Download the antiSMASH database, needed for the prepare-data process
# Prepare the antiSMASH cache and other bits it needs to run properly
# Remove the databases to reduce the container size
RUN download-antismash-databases && antismash --prepare-data && antismash --check-prereqs && rm -rf /usr/local/lib/python3.11/site-packages/antismash/databases

RUN python -c "from antismash.modules.nrps_pks import prepare_data; prepare_data()" && \
    cp /usr/local/lib/python3.11/site-packages/antismash/modules/nrps_pks/data/terminals/*.brawn_cache /usr/local/lib/python3.11/site-packages/antismash/modules/nrps_pks/data/
