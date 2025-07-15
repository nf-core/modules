# Test build locally before making a PR (from project root directory)
#   docker build . -t devcontainer:local -f nf_core/devcontainer/devcontainer.Dockerfile

# Uses mcr.microsoft.com/devcontainers/base:ubuntu
FROM ghcr.io/nextflow-io/training@sha256:97cce091b2c786f8fbd86f470e59d096dff546fe07941cf0e97421b6f95333e2

# Install apptainer via conda
RUN conda install --quiet --yes --update-all --name base \
    conda-forge::apptainer>=1.4.1 && \
    conda clean --all --force-pkgs-dirs --yes

# Change ownership of conda and nextflow bin
RUN chown -R vscode:vscode /usr/local/bin/ && \
    chown -R vscode:vscode /opt/conda

USER vscode
WORKDIR /home/vscode/

# install nf-test to local
RUN wget -qO- https://get.nf-test.com > install_nf_test.sh && \
    bash install_nf_test.sh && \
    mkdir "$HOME/bin" && \
    mv nf-test "$HOME/bin/nf-test" && \
    echo "export PATH=$HOME/bin:$PATH" >> "$HOME/.bashrc"
