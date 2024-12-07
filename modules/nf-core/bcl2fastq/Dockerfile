# syntax=docker/dockerfile:1

# This does everything that we'd want to do to build from the source code
# then we also don't have to download the RPM from Illumina,
# since they're probably not going to make a new version. ;)
FROM zymoresearch/bcl2fastq:2.20

LABEL org.opencontainers.image.authors="Matthias De Smet <matthias.desmet@ugent.be>, Edmund Miller <edmund@nf-co.re>"
LABEL org.opencontainers.image.description="Docker image containing bcl2fastq"
LABEL org.opencontainers.image.version="$BCL2FASTQ_VERSION"
LABEL org.opencontainers.image.documentation="https://github.com/nf-core/modules/blob/master/modules/nf-core/bcl2fastq/README.md"
LABEL org.opencontainers.image.source="https://github.com/nf-core/modules"
LABEL org.opencontainers.image.vendor="nf-core"
LABEL org.opencontainers.image.license="https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-eula.pdf"
# Disclaimer: this container is not provided nor supported by Illumina

# 'ps' command is need by some nextflow executions to collect system stats
# Install procps and clean up
RUN apk --no-cache add \
    bash \
    procps

# Link hostname cmd to fix hardcoded path
RUN ln -s /bin/hostname /usr/bin/hostname

# ENTRYPOINT ["/bin/bash"]
