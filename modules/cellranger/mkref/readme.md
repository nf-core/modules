## Updating the docker container and making a new module release

Cell Ranger is a commercial tool and cannot be distributed. Updating the Cell Ranger version in the container and pushing the update to Dockerhub needs
needs to be done manually.

1. Navigate to the [Cell Ranger download page](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) and download the tar ball of the desired Cell Ranger version with `curl` or `wget`. Place this file in the same folder where the Dockerfile lies.

2. Edit the Dockerfile: update the Cell Ranger version in this line:

    ```bash
    ENV CELLRANGER_VER <VERSION>
    ```

3. Create the container:

    ```bash
    docker build . -t nfcore/cellranger:<VERSION>
    docker push nfcore/cellranger:<VERSION>
    ```
