FROM mambaorg/micromamba:2-amazon2023 AS build
USER root
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.gpu.yml /tmp/environment.gpu.yml
# expose `which` at /usr/bin/which for R (bioconda) post-link scripts; the amazon2023 base image lacks it
RUN CONDA_OVERRIDE_CUDA="99" micromamba install -y -n base -f /tmp/environment.gpu.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba clean -a -y \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END"
RUN micromamba clean --all --yes

FROM ubuntu:24.04 AS prod
ARG MAMBA_ROOT_PREFIX="/opt/conda"
ENV MAMBA_ROOT_PREFIX=$MAMBA_ROOT_PREFIX
COPY --from=build "$MAMBA_ROOT_PREFIX" "$MAMBA_ROOT_PREFIX"
USER root
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
