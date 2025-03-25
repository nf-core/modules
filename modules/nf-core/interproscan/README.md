# Updating the docker container and making a new module release

As of March 2025, the latest version of [InterProScan on Bioconda](https://bioconda.github.io/recipes/interproscan/README.html) is 5.59_91, which is from [October 17, 2022](https://github.com/ebi-pf-team/interproscan/releases/tag/5.59-91.0). This Dockerfile builds a new version of InterProScan.

1. Navigate to the appropriate download page. - [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest): download the tar ball of the desired Cell Ranger version with `curl` or `wget`. Place this file in the same folder where the Dockerfile lies.

2. Edit the Dockerfile. Update the Cell Ranger versions in this line:

```bash
ENV INTERPROSCAN_VER=<VERSION>
```

3. Create and test the container:

You can do `make build` from the Makefile or: 

```bash
docker build . -t quay.io/nf-core/interproscan:<VERSION>
```

4. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

You can do `make push` from the Makefile or:

```bash
docker push quay.io/nf-core/interproscan:<VERSION>
```