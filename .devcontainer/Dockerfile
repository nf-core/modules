# Test build locally before making a PR (from project root directory)
#   docker build . -t devcontainer:local -f nf_core/devcontainer/devcontainer.Dockerfile

# Uses mcr.microsoft.com/devcontainers/base:ubuntu
FROM ghcr.io/nextflow-io/training@sha256:97cce091b2c786f8fbd86f470e59d096dff546fe07941cf0e97421b6f95333e2

# install nf-test to local
RUN wget -qO- https://get.nf-test.com > install_nf_test.sh && \
    bash install_nf_test.sh && \
    mv nf-test "$HOME/.local/bin"

