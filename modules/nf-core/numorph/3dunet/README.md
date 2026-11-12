# Updating the docker container and making a new module release

To update this tool : 2. Edit the Dockerfile in the following lines:

```diff
- ARG NUMORPH_PREPROCESSING_VERSION="1.0.0"
+ ARG NUMORPH_PREPROCESSING_VERSION="<New version>"
```

3. Open a PR against [nf-core/modules](https://github.com/nf-core/modules) directly, push the changes, and the Dockerfile should be built and uploaded to `quay.io/nf-core/numorph-3dunet`!
