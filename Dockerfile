# Test build locally before making a PR (from project root directory)
#   docker build . -t devcontainer:local -f nf_core/devcontainer/devcontainer.Dockerfile

# Uses mcr.microsoft.com/devcontainers/base:ubuntu
FROM ghcr.io/nextflow-io/training@sha256:97cce091b2c786f8fbd86f470e59d096dff546fe07941cf0e97421b6f95333e2

# Add the nf-core source files to the image
COPY --chown=vscode:vscode . /usr/src/nf_core_modules
WORKDIR /usr/src/nf_core_modules

# Install apptainer via conda
RUN conda install --quiet --yes --update-all --name base \
    conda-forge::apptainer>=1.4.1 && \
    conda clean --all --force-pkgs-dirs --yes

# Change ownership of conda and nextflow bin
RUN chown -R vscode:vscode /usr/local/bin/ && \
    chown -R vscode:vscode /opt/conda && \
    chown -R vscode:vscode /workspaces/

USER vscode
