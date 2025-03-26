FROM mambaorg/micromamba:2.0.5

ARG VERSION=5.73-104.0

LABEL software="InterProScan"
LABEL software.version="5.73-104.0"
LABEL description="InterProScan is a search engine for InterPro is ae database which integrates together predictive information about proteins' function."
LABEL website="https://github.com/ebi-pf-team/interproscan"
LABEL license="Apache 2.0"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.lock /tmp/env.lock

RUN micromamba install --name base --yes --file /tmp/env.lock \
    && micromamba clean --all --yes

# This is to install libgomp - installing this with conda doesn't work for bin/cdd/rpsblast
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1==12.2.0-14 \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root

WORKDIR /opt

RUN wget --progress ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/alt/interproscan-core-${VERSION}.tar.gz && \
    wget --progress ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/alt/interproscan-core-${VERSION}.tar.gz.md5 && \
    md5sum -c interproscan-core-${VERSION}.tar.gz.md5 && \
    tar -zxf interproscan-core-${VERSION}.tar.gz && \
    rm interproscan-core-${VERSION}.tar.gz* && \
    mv /opt/interproscan-${VERSION} /opt/interproscan

WORKDIR /opt/interproscan

RUN echo "binary.phobius.pl.path=/opt/interproscan/licensed/phobius/phobius.pl" >> interproscan.properties && \
    echo "binary.signalp.path=/opt/interproscan/licensed/signalp/signalp" >> interproscan.properties && \
    echo "signalp.perl.library.dir=/opt/interproscan/licensed/signalp/lib" >> interproscan.properties && \
    echo "binary.tmhmm.path=/opt/interproscan/licensed/tmhmm/decodeanhmm.Linux_x86_64" >> interproscan.properties && \
    echo "tmhmm.model.path=/opt/interproscan/licensed/tmhmm/TMHMM2.0.model" >> interproscan.properties

ENV PATH="$MAMBA_ROOT_PREFIX/bin:/opt/interproscan/:$PATH"

USER $MAMBA_USER
