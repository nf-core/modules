# Updating the docker container and making a new module release

The GPU build for COMEBin-LLMopt is too large to complete on the standard Seqera container runner, so we have to build a valid container ourselves.

To update the container:

1. Update the `environment.gpu.yaml` to have the latest combin-llmopt version, and compatible versions of CUDA and pytorch-gpu.

2. Create and test the container:

```bash
docker build -f Dockerfile -t quay.io/nf-core/comebin-llmopt:<VERSION> .
```

3. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

```bash
docker push quay.io/nf-core/comebin-llmopt:<VERSION>
```
