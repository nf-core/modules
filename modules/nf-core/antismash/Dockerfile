FROM quay.io/biocontainers/antismash:8.0.0--pyhdfd78af_0

# Download the antiSMASH database, needed for the prepare-data process
# Prepare the antiSMASH cache and other bits it needs to run properly
# Remove the databases to reduce the container size
RUN download-antismash-databases && antismash --prepare-data && antismash --check-prereqs && rm -rf /usr/local/lib/python3.11/site-packages/antismash/databases

RUN python -c "from antismash.modules.nrps_pks import prepare_data; prepare_data()" && \
    cp /usr/local/lib/python3.11/site-packages/antismash/modules/nrps_pks/data/terminals/*.brawn_cache /usr/local/lib/python3.11/site-packages/antismash/modules/nrps_pks/data/
