# Updating the docker container and making a new module release

Universc depends on Cellranger 3.0.2.
Cell Ranger is a commercial tool from 10X Genomics. The container provided for the cellranger nf-core module is not provided nor supported by 10x Genomics. Updating the Cell Ranger versions in the container and pushing the update to Dockerhub needs to be done manually.

1. Navigate to the appropriate download page. - [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/downloads#download-links).
   And copy the wget URL.

```bash
wget -O cellranger-10.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.1.0.tar.gz?Expires=xxxxxx"
```

2. Edit the Dockerfile.
   Update the softwares versions in these lines:

```bash
ENV CELLRANGER_VER=10.1.0
ENV UNIVERSC_VER=1.2.7
ENV MULTIQC_VER=1.35
ENV FASTQ_PAIR_VER=0.3
ENV SAMTOOLS_VER=1.16.1
ENV STAR_VER=2.5.1b
```

and the URL for cellranger

```bash
RUN wget -O cellranger-$CELLRANGER_VER.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-${CELLRANGER_VER}.tar.gz?Expires=XXXX&Key-Pair-Id=XXXXX" \
```

3. Create and test the container:

```bash
wave \
  --containerfile Dockerfile \
  --context . \
  --await
```

4. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

```bash
docker push quay.io/nf-core/universc:<VERSION>
```
