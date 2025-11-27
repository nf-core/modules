#  Use a base image with the required operating system
FROM ubuntu:20.04 AS builder

LABEL authors="Sameesh Kher <sameesh.kher@dkfz-heidelberg.de>" \
    description="Docker image containing Xenium Ranger"

RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y wget \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Download and extract xeniumranger
# The latest version of the xeniumranger tool is downloaded from https://www.10xgenomics.com/support/software/xenium-ranger/downloads
RUN wget -q -O xeniumranger.tar.gz "https://cf.10xgenomics.com/releases/xeniumranger/xeniumranger-4.0.0.tar.gz?Expires=1764227131&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=Be~2wTlK13MnrKJ7d7f1uC5BRr7hE1DzTYEDV7swJCGhz8ey6cNkTat-cbTWO9hI99R6Bp5CMH0aIOok~TEiejXfxC2F-4yJwv1xTEGfP7Oyj87VHQG5p0tfPAXfKK4k3ua5VFdxEjRhPubcsAGXyI6IgIHlHEP1nJTBA2VX4oIE5v7Vw0IYB9pV6-scZs9aKxfOOwI19AoNXP16FYXKqNePRUWe8IHC8NS1ClhoDK5MKjSyYqKtJY0terOKFkmS9kKwancspjp~UKfDPZXK~OiH24dguhWX7PCtrPcSYs9ajiDIYUtoTvGWJmw~p3vvR8wk6Am03W0d95yTmMHqeg__" && \
    tar -xzvf xeniumranger.tar.gz && \
    rm xeniumranger.tar.gz

# multistage to reduce image size
FROM ubuntu:20.04

# Set environment variables
ENV PATH="/xeniumranger-xenium4.0/bin:$PATH"

# copy over xenium from builder
COPY --from=builder /xeniumranger-xenium4.0 /xeniumranger-xenium4.0/
