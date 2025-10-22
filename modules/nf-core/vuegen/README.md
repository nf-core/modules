# VueGen module

## Images using wave

Cached images for 0.5.1 version available on wave, see 
[docs](https://docs.seqera.io/wave#seqera-containers---the-community-container-registry)

platform    | version     | image
--------    | -------     | -----
docker      | linux/amd64 | community.wave.seqera.io/library/vuegen_python:236414fc5cfce774
singularity | linux/amd64 | oras://community.wave.seqera.io/library/vuegen_python:1adb57ecbfa02088
docker      | linux/arm64 | community.wave.seqera.io/library/vuegen_python:a1dc6490f5045706
singularity | linux/arm64 | oras://community.wave.seqera.io/library/vuegen_python:64670fcadedf151f

## Custom docker file

- check Docker best practices with hadolint, see
 [docker-docs](https://docs.docker.com/build/building/best-practices/)


```bash
$ hadolint Dockerfile
Dockerfile:7 DL3008 warning: Pin versions in apt get install. Instead of `apt-get install <package>` use `apt-get install <package>=<version>`
Dockerfile:24 DL3059 info: Multiple consecutive `RUN` instructions. Consider consolidation.
Dockerfile:31 SC2211 warning: This is a glob used as a command name. Was it supposed to be in ${..}, array, or is it missing quoting?
Dockerfile:31 DL4006 warning: Set the SHELL option -o pipefail before RUN with a pipe in it. If you are using /bin/sh in an alpine image or if your shell is symlinked to busybox then consider explicitly setting your SHELL to /bin/ash, or disable this check
```
