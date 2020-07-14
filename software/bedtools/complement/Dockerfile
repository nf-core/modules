FROM nfcore/base
LABEL authors="Jose Espinosa-Carrasco" \
    description="Docker image containing all requirements for nf-core/modules/bedtools/complementbed"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-bedtools-complementbed/bin:$PATH
