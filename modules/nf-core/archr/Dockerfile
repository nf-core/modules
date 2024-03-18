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
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


# install R packages
RUN R -q -e "install.packages('BiocManager')"
RUN R -q -e "install.packages('devtools')"
RUN R -q -e "library(devtools)"
RUN R -q -e "devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())"
RUN R -q -e "library(ArchR)"
RUN R -q -e "ArchR::installExtraPackages()"

# maybe install hg19/hg28,mm9,mm10 BSgenomes
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")

# this might not give error if installation fails

