# Updating the docker container and making a new module release

Cell Ranger ATAC is a commercial tool from 10X Genomics. The container provided for the cellranger nf-core module is
not provided nor supported by 10x Genomics. Updating the Cell Ranger ATAC versions in the container and pushing the
update to Dockerhub needs to be done manually.

1. Navigate to the appropriate download page. - [Cell Ranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation): download the tar ball of the desired Cell Ranger ATAC version with `curl` or `wget`. Place this
   file in the same folder where the Dockerfile lies.

2. Edit the Dockerfile. Update the Cell Ranger versions in this line:

   ```bash
   ENV CELLRANGERATAC_VER=<VERSION>
   ```

3. Create and test the container:

   ```bash
   docker build . -t nfcore/cellrangeratac:<VERSION>
   ```

4. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member
   to do so.

   ```bash
   docker push nfcore/cellrangeratac:<VERSION>
   ```
