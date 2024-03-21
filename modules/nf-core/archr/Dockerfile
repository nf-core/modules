# dockerfile with R, biocManager, and archer
FROM rocker/r-ubuntu:20.04


# Override labels from base image
LABEL authors="xxx" \
    description="Docker image containing archR"

# install apt-get packages
RUN apt-get -y update && apt-get install --no-install-recommends -y \
    build-essential \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libzmq3-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libgsl-dev \ 
    && \ 
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


# install R packages
RUN R -q -e "install.packages('BiocManager')"
RUN R -q -e "install.packages('remotes')"
RUN R -q -e "install.packages('devtools')"
RUN R -q -e "library(devtools)"
#RUN installGithub.r GreenleafLab/ArchR

RUN R -q -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"
RUN R -q -e "library(ArchR)"
RUN R -q -e "ArchR::installExtraPackages()"

# support for different genomes
RUN R -e "BiocManager::install(c('BiocVersion', 'BSgenome.Hsapiens.UCSC.hg38', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Mmusculus.UCSC.mm9', 'BSgenome.Mmusculus.UCSC.mm10'), update=TRUE, ask=FALSE)"


# this might not give error if installation fails

