# Updating the docker container and making a new module release

The CellPheno NIS container is built from the [`Dockerfile`](./Dockerfile) in this
directory. It clones a pinned release of the source repository, compiles the
C++/LibTorch/CUDA executable, and installs it on `PATH` as `cellpheno-nis`. The module
reports the tool version at runtime via `cellpheno-nis --version`.

To update this tool:

1. In the [`Dockerfile`](./Dockerfile), pin the source to the new release tag:

   ```diff
   - ARG NIS_SOURCE_REF=v1.0.3
   + ARG NIS_SOURCE_REF=<new release tag>
   ```

2. In [`main.nf`](./main.nf), bump the `container` tag to match:

   ```diff
   - container "quay.io/nf-core/cellpheno-nis:1.0.3"
   + container "quay.io/nf-core/cellpheno-nis:<new version>"
   ```

3. Open a PR against [nf-core/modules](https://github.com/nf-core/modules) and push the
   changes; the `Dockerfile` is built and uploaded to `quay.io/nf-core/cellpheno-nis`.
