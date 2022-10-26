# Updating the docker container and making a new module release

Cell Ranger is a commercial tool from 10X Genomics. The container provided for the cellranger nf-core module is not provided nor supported by 10x Genomics. Updating the Cell Ranger versions in the container and pushing the update to Dockerhub needs to be done manually.

1. Navigate to the appropriate download page. - [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest): download the tar ball of the desired Cell Ranger version with `curl` or `wget`. Place this file in the same folder where the Dockerfile lies.

2. Edit the Dockerfile. Update the Cell Ranger versions in this line:

   ```bash
   ENV CELLRANGER_VER=<VERSION>
   ```

3. Create and test the container:

   ```bash
   docker build . -t nfcore/cellranger:<VERSION>
   ```

4. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

   ```bash
   docker push nfcore/cellranger:<VERSION>
   ```
