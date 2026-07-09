# Container image for the CellPheno NIS executable (modules/nf-core/cellpheno/nis).
#
# Self-contained build: clones the source repository, compiles the C++/LibTorch/OpenCV/
# CUDA binary and ships it on PATH as `main`. nf-core builds this Dockerfile from the
# module directory and uploads it to quay.io/nf-core/cellpheno-nis (see README.md), which
# the module's `container` line references.
#
# GPU container: run with `--gpus all` (docker) / `--nv` (singularity).

# CUDA 12.1 (not 12.8) is chosen for driver compatibility: a 12.1 image runs on any
# NVIDIA driver >= 530 (vs. >= ~570 for 12.8), covering far more clusters/CI runners while
# still supporting every GPU this tool targets. The binary has no custom .cu kernels
# (CMakeLists compiles only .cpp), so GPU-arch coverage comes entirely from LibTorch's
# prebuilt fat binary: cu121 spans sm_50..sm_90 (Maxwell->Hopper, incl. Ada sm_89). This
# also matches the cellpheno-postproc container (torch 2.5.1+cu121).
ARG CUDA_VERSION=12.1.1
ARG UBUNTU_VERSION=22.04

# ---------------------------------------------------------------------------
# Stage 1: build the NIS binary from source
# ---------------------------------------------------------------------------
# Base without cuDNN: LibTorch bundles its own cuDNN in libtorch/lib (first on
# LD_LIBRARY_PATH), so we avoid pinning a base cuDNN tag.
FROM nvidia/cuda:${CUDA_VERSION}-devel-ubuntu${UBUNTU_VERSION} AS build

# LibTorch CUDA build matching the base CUDA version.
ARG LIBTORCH_URL=https://download.pytorch.org/libtorch/cu121/libtorch-cxx11-abi-shared-with-deps-2.5.1%2Bcu121.zip
# Source of the NIS C++ executable, pinned to a tagged release for reproducibility
# (not a moving branch). v1.0.0 == commit 5420a1cc4b1a7dc8faa71fd7cd36d271fb2f8364.
ARG NIS_SOURCE_REPO=https://github.com/Chrisa142857/Lightsheet_microscopy_image_3D_nuclei_instance_segmentation
ARG NIS_SOURCE_REF=v1.0.0

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        wget \
        unzip \
        ca-certificates \
        libopencv-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN wget -q -O libtorch.zip "${LIBTORCH_URL}" \
    && unzip -q libtorch.zip \
    && rm libtorch.zip

RUN git clone --depth 1 --branch "${NIS_SOURCE_REF}" "${NIS_SOURCE_REPO}" /src
WORKDIR /src/cpp/build
RUN cmake -DCMAKE_PREFIX_PATH=/opt/libtorch -DCMAKE_BUILD_TYPE=Release .. \
    && cmake --build . --config Release -j "$(nproc)"

# ---------------------------------------------------------------------------
# Stage 2: slim runtime image
# ---------------------------------------------------------------------------
FROM nvidia/cuda:${CUDA_VERSION}-runtime-ubuntu${UBUNTU_VERSION} AS runtime

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        libopencv-core4.5d \
        libopencv-imgcodecs4.5d \
        libopencv-imgproc4.5d \
        libgomp1 \
        procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/libtorch /opt/libtorch
COPY --from=build /src/cpp/build/main /usr/local/bin/main

# Single authoritative version source. The module reports the NIS version via
# `eval("cat /usr/local/share/cellpheno-nis/VERSION")`, so the container tag and the
# reported version stay in lockstep. Bump together with the `container` tag in main.nf.
ARG CELLPHENO_NIS_VERSION=1.0.0
RUN mkdir -p /usr/local/share/cellpheno-nis \
    && echo "${CELLPHENO_NIS_VERSION}" > /usr/local/share/cellpheno-nis/VERSION

ENV LD_LIBRARY_PATH=/opt/libtorch/lib:${LD_LIBRARY_PATH}

LABEL org.opencontainers.image.authors="Ziquan Wei (@Chrisa142857)" \
      org.opencontainers.image.title="cellpheno-nis" \
      org.opencontainers.image.description="Whole-brain 3D nuclei instance segmentation for lightsheet microscopy" \
      org.opencontainers.image.source="https://github.com/Chrisa142857/Lightsheet_microscopy_image_3D_nuclei_instance_segmentation" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="1.0.0"

ENTRYPOINT []
CMD ["main", "--help"]
