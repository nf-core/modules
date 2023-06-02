# Updating the docker container and making a new module release

Space Ranger is a commercial tool by 10X Genomics. The container provided for the spaceranger nf-core module is not provided nor supported by 10x Genomics. Updating the Space Ranger version in the container and pushing the update to Dockerhub needs to be done manually.

1. Navigate to the [Space Ranger download page](https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest) and download the tar ball of the desired Space Ranger version with `curl` or `wget`. Place this file in the same folder where the Dockerfile lies.

2. Edit the Dockerfile: update the Cell Ranger version in this line:

   ```bash
   ENV SPACERANGER_VER <VERSION>
   ```

3. Create the container:

   ```bash
   docker build . -t nfcore/spaceranger:<VERSION>
   docker push nfcore/spaceranger:<VERSION>
   ```
