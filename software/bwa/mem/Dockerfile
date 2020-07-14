FROM nfcore/base
LABEL authors="Jeremy Guntoro" \
    description="Docker image containing all requirements for nf-core/modules/bwa/mem module"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-bwa-mem/bin:$PATH
