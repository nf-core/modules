# syntax=docker/dockerfile:1

FROM continuumio/miniconda3:4.9.2
LABEL authors="Gisela Gabernet <gisela.gabernet@gmail.com>,Edmund Miller <edmund.a.miller@gmail.com>" \
    description="Docker image containing Space Ranger"
# Disclaimer: this container is not provided nor supported by 10x Genomics.

# Install procps and curl and clean apt cache
RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y procps curl \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Copy pre-downloaded spaceranger file
ENV SPACERANGER_VER 3.1.2
ADD "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-3.1.2.tar.gz?Expires=1732608367&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=S4jcGCs5H0wLFfREhSc~AfuPIOXE~PW96kX9o2EyxZnmn1goUecgfEWTj67fL1DnZJBIa06kFdUdNpsMn2ustiSWOfXEgjlYQnxIcNnQIiEBGoQTHWphvx3LEQ6wtZnkWS80P6IcE0HJkIsgy04t6Sohih5cxY4jgytYsrAfZDYr5G3KKFwTfCKmhzMaXqW635yPbyQ8xEcQHK0QwviAx8-EFq-PE8UzC4QgUKi2MW-ivcfZkSDSfF8C3s7SgwDXIGIWv52mzeszenxMjN4KrWQotZ7ZpktzI0Vfpz0dNC17dQeDQUHj4LuNYbdh3RqsPKtqu3wjCe2Q7KiyoWnmaw__" \
    /opt/spaceranger-$SPACERANGER_VER \
    --checksum 949088b1e2cc5bd03079b60a65b41d39

# Install spaceranger
RUN \
    cd /opt && \
    tar -xzvf spaceranger-$SPACERANGER_VER.tar.gz && \
    export PATH=/opt/spaceranger-$SPACERANGER_VER:$PATH && \
    ln -s /opt/spaceranger-$SPACERANGER_VER/spaceranger /usr/bin/spaceranger && \
    rm -rf /opt/spaceranger-$SPACERANGER_VER.tar
