FROM nfcore/base
LABEL authors="Niclas Kildegaard Nielsen and Veit Schwämmle" \
      description="Docker image containing all requirements for maxquant pipeline"

COPY environment.yml /
# Create the environment:
RUN conda env create -f /environment.yml && conda clean -a
# Make RUN commands use the new environment:

# Activates the env for the program to run in
RUN echo "source activate nf-core-maxquant" > ~/.bashrc
ENV PATH /opt/conda/envs/nf-core-maxquant/bin:$PATH

# Changes the Work dir to the SDRF file area
WORKDIR /
# ENTRYPOINT ["/bin/bash"]
