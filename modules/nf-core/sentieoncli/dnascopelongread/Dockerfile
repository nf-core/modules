FROM debian:stable-20260112-slim AS downloader

ARG SENTIEON_VERSION=202503.02
ARG SENTIEON_CLI_VERSION=1.5.2

RUN test -n "$SENTIEON_VERSION"

LABEL container.base.image="debian:stable-20260112-slim" \
    software.version="${SENTIEON_VERSION}" \
    software.website="https://www.sentieon.com/"

# Download the software from the permalink
RUN apt-get update && apt-get install -y curl
RUN mkdir -p /opt/sentieon/
RUN curl -L "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${SENTIEON_VERSION}.tar.gz" | tar -zxf - -C /opt/sentieon/

# Build gnuplot
RUN apt-get update && apt-get install -y curl make gcc gzip libgd-dev libcairo2-dev libpango1.0-dev && \
    curl -L "https://downloads.sourceforge.net/project/gnuplot/gnuplot/6.0.2/gnuplot-6.0.2.tar.gz" | \
    tar -zxf - && \
    cd gnuplot-6.0.2 && \
    ./configure && \
    make install

# Build the container
FROM debian:stable-20260112-slim
ENV SENTIEON_VERSION=$SENTIEON_VERSION

COPY --from=downloader /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION} /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION}
COPY --from=downloader /usr/local/bin/gnuplot /usr/local/bin/gnuplot
COPY --from=downloader /usr/local/libexec/gnuplot/6.0/gnuplot_x11 /usr/local/libexec/gnuplot/6.0/gnuplot_x11
COPY --from=downloader /usr/local/share/gnuplot/6.0 /usr/local/share/gnuplot/6.0
CMD ["/bin/bash"]

# Install jemalloc as the recommended memory allocation tool, see https://support.sentieon.com/appnotes/jemalloc/
# Install procps for process monitoring
RUN apt-get update && apt-get install -y \
    libjemalloc2 \
    procps \
    libgd-dev \
    libcairo2-dev \
    libpango1.0-dev \
    locales \
    git \
    curl \
    ca-certificates \
    bzip2 \
    isal \
    && sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen \
    && dpkg-reconfigure --frontend=noninteractive locales \
    && echo "LANG=en_US.UTF-8" | tee /etc/default/locale \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives

ENV LC_ALL=en_US.UTF-8
ENV SENTIEON_INSTALL_DIR=/opt/sentieon/sentieon-genomics-$SENTIEON_VERSION
ENV PATH=${SENTIEON_INSTALL_DIR}/bin:${PATH}
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2

ARG MICROMAMBA_VERSION=1.5.8-0
ARG MICROMAMBA_ARCH=linux-64
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=${MAMBA_ROOT_PREFIX}/bin:${PATH}

RUN set -eux; \
    mkdir -p "$MAMBA_ROOT_PREFIX"; \
    curl -fsSL "https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-${MICROMAMBA_ARCH}" -o /usr/local/bin/micromamba; \
    chmod 755 /usr/local/bin/micromamba

ARG SENTIEON_CLI_REPO=https://github.com/Sentieon/sentieon-cli.git

RUN set -eux; \
    micromamba create -y -n sentieon-cli -c conda-forge -c bioconda \
    python=3.11 \
    pip=24.0 \
    samtools=1.22 \
    bcftools=1.22 \
    htslib=1.22 \
    pysam=0.23.0 \
    click=8.1.7 \
    pyyaml=6.0.1 \
    pandas=2.2.2 \
    requests=2.32.3 \
    bedtools=2.31.1 \
    mosdepth=0.3.12 \
    hificnv=1.0.1 \
    packaging=24.0 \
    colorlog=6.9.0 \
    importlib_resources=6.4.5; \
    micromamba clean --all --yes

RUN set -eux; \
    micromamba run -n sentieon-cli pip install --no-deps "git+${SENTIEON_CLI_REPO}@${SENTIEON_CLI_VERSION}";

ENV SENTIEON_CLI_ENV=${MAMBA_ROOT_PREFIX}/envs/sentieon-cli
ENV PATH=${SENTIEON_CLI_ENV}/bin:${MAMBA_ROOT_PREFIX}/bin:${SENTIEON_INSTALL_DIR}/bin:${PATH}
ENV CONDA_DEFAULT_ENV=sentieon-cli

# A default jemalloc configuration that should work well for most use-cases, see http://jemalloc.net/jemalloc.3.html
ENV MALLOC_CONF=metadata_thp:auto,background_thread:true,dirty_decay_ms:30000,muzzy_decay_ms:30000
