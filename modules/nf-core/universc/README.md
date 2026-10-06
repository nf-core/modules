# Updating the docker container and making a new module release

Universc depends on Cellranger 3.0.2.
This version is already available through tomkellygenetics/cellranger_clean:3.0.2.9002.

1. Edit the Dockerfile. Update the Cell Ranger and universc versions:

```bash
ENV CELLRANGER_VER=<VERSION>
```

2. Create and test the container:

```bash
wave \
  --containerfile Dockerfile \
  --context . \
  --await
```

3. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

```bash
docker push quay.io/nf-core/universc:<VERSION>
```
