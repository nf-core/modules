# Updating the docker container and making a new module release

The CellPheno NIS container is built from the [`Dockerfile`](./Dockerfile) in this
directory. It clones a pinned release of the source repository, compiles the
C++/LibTorch/CUDA executable, and installs it on `PATH` as `cellpheno-nis`.

To update this tool:

1. In the [`Dockerfile`](./Dockerfile), pin the source to the new release tag and bump the
   reported version:

   ```diff
   - ARG NIS_SOURCE_REF=v1.0.2
   + ARG NIS_SOURCE_REF=<new release tag>
   - ARG CELLPHENO_NIS_VERSION=1.0.2
   + ARG CELLPHENO_NIS_VERSION=<new version>
   ```

2. In [`main.nf`](./main.nf), bump the `container` tag and the version string in the
   `versions` topic to match:

   ```diff
   - container "quay.io/nf-core/cellpheno-nis:1.0.2"
   + container "quay.io/nf-core/cellpheno-nis:<new version>"
   ...
   - tuple val("${task.process}"), val('cellpheno-nis'), val('1.0.2'), topic: versions, ...
   + tuple val("${task.process}"), val('cellpheno-nis'), val('<new version>'), topic: versions, ...
   ```

3. Open a PR against [nf-core/modules](https://github.com/nf-core/modules) and push the
   changes; the `Dockerfile` is built and uploaded to `quay.io/nf-core/cellpheno-nis`.
