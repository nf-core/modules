# Updating the docker container and making a new module release

Xenium Ranger is a commercial tool from 10X Genomics. The container provided for the xeniumranger nf-core module is not provided nor supported by 10x Genomics. Updating the Xenium Ranger versions in the container and pushing the update to Dockerhub needs to be done manually.

1. Navigate to the appropriate download page. - [Xenium Ranger](https://www.10xgenomics.com/support/software/xenium-ranger/downloads): download the tar ball of the desired Xenium Ranger version with `curl` or `wget`. Place this file in the same folder where the Dockerfile lies.

2. Edit the Dockerfile. Update the Xenium Ranger versions in this line:
   (current version for xenium ranger is fixated at [3.0.1](https://www.10xgenomics.com/support/software/xenium-ranger/downloads) in the dockerfile)

3. Create and test the container:

```bash
docker build . -t quay.io/nf-core/xeniumranger:<VERSION>
```

4. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

```bash
docker push quay.io/nf-core/xeniumranger:<VERSION>
```
