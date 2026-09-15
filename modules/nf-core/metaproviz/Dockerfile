# Shared container for the MetaProViz module family:
# metaproviz/poolestimation, metaproviz/processing, metaproviz/internalstandard
# (they all use the same R + MetaProViz stack, so one image serves all three).
#
# MetaProViz is a Bioconductor package (release 4.0.0, BioC 3.23) but
# not on bioconda. So we still bake it into this custom image, but install
# it from the offical Bioconductor release via BiocManager rather than
# pinning an arbitrary GitHub commit — more reproducible, citable, versioned.
#
# BioC 3.23 is paired specifically with R 4.6.x.
#
# No entrypoint (nf-core injects its own command). rocker/r-ver:4.6.1 ships
# both linux/amd64 and linux/arm64 as its base — but that alone doesn't make
# YOUR built image multi-arch (see build command below).
#
# Build once + test + push (Docker Desktop works fine locally, incl. on
# macOS — no remote/Ona environment needed for this step):
#   cd modules/nf-core/metaproviz
#
#   # quick local test build (single-arch, whatever machine you're on):
#   docker build -t metaproviz-test .
#   docker run --rm metaproviz-test pool_estimation.R --help || true
#
#   # real build + push — use buildx --platform, NOT plain `docker build`.
#   # A plain build only builds for whichever machine runs it (e.g.
#   # arm64-only if built on an Apple Silicon Mac), which would break for
#   # anyone on amd64 (e.g. Gernot on Windows, where Docker Desktop runs
#   # over WSL2 = amd64). --platform covers both in one pushed tag:
#   docker buildx build --platform linux/amd64,linux/arm64 \
#       -t ghcr.io/saezlab/metaproviz:1.0.0 --push .

FROM rocker/r-ver:4.6.1

# zlib1g-dev, libglpk-dev, libuv1-dev: required transitively by MetaProViz's
# Bioconductor dependency chain (XVector needs zlib.h to compile; igraph's
# compiled .so needs libglpk at runtime for CARNIVAL; fs's compiled .so needs
# libuv at runtime for OmnipathR). Without these, BiocManager::install()
# below reports errors per-package but still exits 0 overall, so a broken
# image builds "successfully" 
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        zlib1g-dev libglpk-dev libuv1-dev \
    && rm -rf /var/lib/apt/lists/*

# options(timeout=600): R's default download timeout is 60s. A multi-arch
# build failed here with "download from '.../PACKAGES' failed" — traced it
# to the actual network, not Docker/emulation: a plain curl to the same
# mirror only received ~70KB in a full 60s (a PACKAGES index is several MB),
# so R's default timeout was simply too short for this network's current
# throughput, on either CRAN mirror tried. 600s gives slow-but-working
# downloads enough room without masking a genuinely broken connection.
RUN R -e "options(timeout = 600); install.packages(c('remotes','dplyr','tidyr','tibble','magrittr','ggplot2','rlang','ggfortify','base64enc','htmltools'), repos='https://cloud.r-project.org')" \
 && R -e "options(timeout = 600); install.packages('BiocManager'); BiocManager::install(version = '3.23', update = FALSE, ask = FALSE)" \
 && R -e "options(timeout = 600); BiocManager::install(c('SummarizedExperiment', 'S4Vectors', 'MetaProViz'), update = FALSE, ask = FALSE)"

# Guide's cache-folder rule: tool caches/logs must land in a writable dir, not a
# read-only $HOME. MetaProViz / OmnipathR write log + cache dirs (metaproviz-log/,
# omnipathr-log/) and fetch reference tables on first use.
ENV HOME=/tmp
ENV R_USER_CACHE_DIR=/tmp/rcache
ENV XDG_CACHE_HOME=/tmp/xdg-cache
RUN mkdir -p /tmp/rcache /tmp/xdg-cache && chmod -R 777 /tmp
