# Updating the Docker image

So far no public docker image is available for MuSE. This might change in the future, see Issue [#21](https://github.com/wwylab/MuSE/issues/21) on the official repository. Until then we have to build it ourselves using the provided Dockerfile.

1. Build docker image using provided Dockerfile.

   ```bash
   docker build -t muse:<VERSION> --platform linux/amd64 .
   ```

2. Access rights are needed to push the container to the Dockerhub nf-core organization, please ask a core team member to do so.

   ```bash
   docker tag muse:<VERSION> quay.io/nf-core/muse:<VERSION>
   docker push quay.io/nf-core/muse:<VERSION>
   ```

3. Make the image public.
