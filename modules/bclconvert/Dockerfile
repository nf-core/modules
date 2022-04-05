# Dockerfile to create container with bcl-convert
# Push to nfcore/bclconvert:<VER>

FROM debian:bullseye-slim
LABEL authors="Matthias De Smet <matthias.desmet@ugent.be>" \
    description="Docker image containing bcl-convert"
# Disclaimer: this container is not provided nor supported by Illumina

COPY bcl-convert /usr/local/bin/bcl-convert
RUN chmod +x /usr/local/bin/bcl-convert
