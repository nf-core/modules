FROM rocker/r-ver:4.6.1

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        zlib1g-dev libglpk-dev libuv1-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "options(timeout = 600); install.packages(c('remotes','dplyr','tidyr','tibble','magrittr','ggplot2','rlang','ggfortify','base64enc','htmltools'), repos='https://cloud.r-project.org')" \
 && R -e "options(timeout = 600); install.packages('BiocManager'); BiocManager::install(version = '3.23', update = FALSE, ask = FALSE)" \
 && R -e "options(timeout = 600); BiocManager::install(c('SummarizedExperiment', 'S4Vectors', 'MetaProViz'), update = FALSE, ask = FALSE)"

ENV HOME=/tmp
ENV R_USER_CACHE_DIR=/tmp/rcache
ENV XDG_CACHE_HOME=/tmp/xdg-cache
RUN mkdir -p /tmp/rcache /tmp/xdg-cache && chmod -R 777 /tmp
