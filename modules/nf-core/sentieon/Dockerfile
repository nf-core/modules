FROM debian:bullseye-20220328-slim AS downloader

ARG SENTIEON_VERSION
RUN test -n "$SENTIEON_VERSION"

LABEL container.base.image="debian:bullseye-20220328-slim" \
    software.version="${SENTIEON_VERSION}" \
    software.website="https://www.sentieon.com/"

# Download the software from the permalink
RUN apt-get update && apt-get install -y curl && \
    mkdir -p /opt/sentieon/ && \
    curl -L "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${SENTIEON_VERSION}.tar.gz" | \
    tar -zxf - -C /opt/sentieon/

# Build the container
FROM debian:bullseye-20220328-slim
ARG SENTIEON_VERSION
ENV SENTIEON_VERSION=$SENTIEON_VERSION

COPY --from=downloader /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION} /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION}
CMD ["/bin/bash"]

# Install jemalloc as the recommended memory allocation tool, see https://support.sentieon.com/appnotes/jemalloc/
# Install procps for process monitoring
RUN apt-get update && apt-get install -y \
    libjemalloc2 \
    procps \
    && rm -rf /var/lib/apt/lists/*

ENV SENTIEON_INSTALL_DIR=/opt/sentieon/sentieon-genomics-$SENTIEON_VERSION
ENV PATH $SENTIEON_INSTALL_DIR/bin/:$PATH
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2

# A default jemalloc configuration that should work well for most use-cases, see http://jemalloc.net/jemalloc.3.html
ENV MALLOC_CONF=metadata_thp:auto,background_thread:true,dirty_decay_ms:30000,muzzy_decay_ms:30000
