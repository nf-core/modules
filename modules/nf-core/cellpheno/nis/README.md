# Updating the docker container and making a new module release

The CellPheno NIS container is built from the [`Dockerfile`](./Dockerfile) in this
directory. It clones the source repository and compiles the C++/LibTorch/CUDA executable,
then ships it on `PATH` as `main`.

To update this tool:

1. Edit the version (and, for a fixed build, the pinned source ref) in the `Dockerfile`:

   ```diff
   - ARG CELLPHENO_NIS_VERSION=1.0.0
   + ARG CELLPHENO_NIS_VERSION=<new version>
   ```

   and bump the `container` tag in [`main.nf`](./main.nf) to match.

2. Open a PR against [nf-core/modules](https://github.com/nf-core/modules), push the
   changes, and the `Dockerfile` will be built and uploaded to
   `quay.io/nf-core/cellpheno-nis`.
