# Updating the docker container and making a new module release

dorado is distributed by Oxford Nanopore Technologies under the Oxford Nanopore Technologies PLC. Public License Version 1.0 and is not available on Bioconda. The container provided for the dorado nf-core modules is not provided nor supported by Oxford Nanopore Technologies. Updating the dorado version in the container and pushing the update to quay.io needs to be done manually.

1. Find the release on the [dorado releases page](https://github.com/nanoporetech/dorado/releases) and compute the checksum of its `linux-x64` tarball:

```bash
curl -fsSL https://cdn.oxfordnanoportal.com/software/analysis/dorado-<VERSION>-linux-x64.tar.gz | sha256sum
```

2. Edit the Dockerfile. Update the version and the checksum in these lines:

```bash
ARG DORADO_VER=<VERSION>
ARG DORADO_SHA256=<SHA256>
```

3. Create and test the container. The build downloads the tarball itself:

```bash
docker build . -t quay.io/nf-core/dorado:<VERSION>
docker run --rm quay.io/nf-core/dorado:<VERSION> dorado --version
```

4. Access rights are needed to push the container to the quay.io nf-core organisation, please ask a core team member to do so.

```bash
docker push quay.io/nf-core/dorado:<VERSION>
```
