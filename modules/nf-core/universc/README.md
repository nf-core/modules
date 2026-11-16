# Updating the docker container and making a new module release

Cell Ranger is a commercial tool from 10X Genomics. The container provided for the universc nf-core module bundles Cell Ranger and is not provided, endorsed, or supported by 10x Genomics.
**Do not redistribute this container or its contents unless you are explicitly authorized to do so under the applicable 10x Genomics license or agreement**. Any authorization applicable to the nf-core organization is limited to the scope of that agreement and must not be assumed to permit independent redistribution or publication elsewhere. Before building or publishing an updated image, confirm that the intended distribution is covered by the applicable authorization.
Updating the Cell Ranger or UniverSC versions in the container and publishing the update needs to be done manually.
UniverSC and Cell Ranger are separate software components with their own licensing terms. Ensure that the applicable terms for both components are respected.

1. Navigate to the appropriate download page. - [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/downloads): download the tar ball of the desired Cell Ranger version with `curl` or `wget`. Place this file in the same folder where the Dockerfile lies.

2. Edit the Dockerfile. Update the versions in these lines:

```bash
ENV CELLRANGER_VER=<VERSION>
ENV UNIVERSC_VER=<VERSION>
```

`universc.patch` adapts UniverSC to the bundled Cell Ranger version. Check it still applies, and update or drop it, when changing either version.

3. Create and test the container:

```bash
docker build -f Dockerfile -t quay.io/nf-core/universc:<VERSION> .
```

4. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so. Before pushing, confirm that publishing the image containing Cell Ranger is covered by the applicable authorization.

```bash
docker push quay.io/nf-core/universc:<VERSION>
```
