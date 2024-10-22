# Updating the docker container and making a new module release

Bcl2fastq2 and Cell Ranger Arc are commercial tools from Illumina and 10X Genomics, respectively. The container provided for the cellranger nf-core module is not provided nor supported by either Illumina or 10x Genomics. Updating the bcl2fastq2 or Cell Ranger Arc versions in the container and pushing the update to Dockerhub needs to be done manually.

1. Navigate to the appropriate download pages. - [bcl2fastq2](https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html): download the linux rpm installer of the desired bcl2fastq2 version with `curl` or `wget`. Place this file in the same folder where the Dockerfile lies. - [Cell Ranger Arc](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/installation): download the tar ball of the desired Cell Ranger Arc version with `curl` or `wget`. Place this file in the same folder where the Dockerfile lies.

2. Edit the Dockerfile. Update the bcl2fastq2 and Cell Ranger Arc versions in this line:

```bash
ENV BCL2FASTQ2_VER=<VERSION> \
    ARC_VER=<VERSION>
```

3. Create and test the container:

```bash
docker build . -t quay.io/nf-core/cellranger-arc-mkfastq:<CELLRANGERARC_VERSION>
```

4. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

```bash
docker push quay.io/nf-core/cellranger-arc-mkfastq:<CELLRANGERARC_VERSION>
```
