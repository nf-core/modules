# Updating the docker container and making a new module release

To update this tool :

1. Navigate tho the Numorph-toolkit [repository](https://github.com/qbic-pipelines/Numorph-toolkit/releases), go to the latest release and get the link and sha256 checksum.
2. Edit the Dockerfile in the following lines:

```diff
- ARG NUMORPH_PREPROCESSING_VERSION="1.0.0"
+ ARG NUMORPH_PREPROCESSING_VERSION="<New version>"
- ARG NUMORPH_PREPROCESSING_URL="https://github.com/qbic-pipelines/Numorph-toolkit/releases/download/1.0.0/numorph_preprocessing"
+ ARG NUMORPH_PREPROCESSING_UR="<New release>"
- ARG NUMORPH_PREPROCESSING_SHA256="2adecf1af4a7d362eca7ab2876e017034c4b155df10fa3adcca409a4ac47c976"
+ ARG NUMORPH_PREPROCESSING_SHA256="<Manual SHA256 of binary file>"
```

3. Open a PR against [nf-core/modules](https://github.com/nf-core/modules) directly, push the changes, and the Dockerfile should be built and uploaded to `quay.io/nf-core/numorph_preprocessing`!
