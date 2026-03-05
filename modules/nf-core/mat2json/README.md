# Updating the docker container and making a new module release
To update this tool :
1. Navigate tho the mat2json [repository](https://github.com/qbic-pipelines/mat2json), go to the latest release and get the link and sha256 checksum.
2. Edit the Dockerfile in the following lines:
```diff
- ARG MAT2JSON_VERSION="1.0.0"
+ ARG MAT2JSON_VERSION="<New version>"
- ARG MAT2JSON_URL="https://github.com/qbic-pipelines/mat2json/releases/download/v1.0.0/mat2json"
+ ARG MAT2JSON_UR="<New release>"
- ARG MAT2JSON_SHA256="216cac1716c761ebcdde415726762ec462cc79eda7daf3e320549b8765f097e5"
+ ARG MAT2JSON_SHA256="<Manual SHA256 of binary file>"
```
3. Open a PR against [nf-core/modules](https://github.com/nf-core/modules) directly, push the changes, and the Dockerfile should be built and uploaded to `quay.io/nf-core/mat2json`!