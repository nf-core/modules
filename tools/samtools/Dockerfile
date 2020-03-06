FROM nfcore/base:1.7
LABEL authors="phil.ewels@scilifelab.se" \
      description="Docker image for nf-core modules samtools"

# foobar
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-modules-samtools/bin:$PATH
