# UniverSC

## Single-cell processing across technologies

UniverSC is an open-source single-cell pipeline that runs across platforms on various technologies.

## Maintainers

Tom Kelly (RIKEN, IMS)

Kai Battenberg (RIKEN CSRS/IMS)

Contact: <first name>.<family name>[at]riken.jp

## Implementation

This container runs Cell Ranger v3.0.2 installed from source on MIT License on GitHub with
modifications for compatibility with updated dependencies. All software is installed from
open-source repositories and available for reuse.

It is _not_ subject to the 10X Genomics End User License Agreement (EULA).
This version allows running Cell Ranger v3.0.2 on data generated from any experimental platform
without restrictions. However, updating to newer versions on Cell Ranger subject to the
10X EULA is not possible without the agreement of 10X Genomics.

To comply with licensing and respect 10X Genomics Trademarks, the 10X Genomics logo
has been removed from HTML reports, the tool has been renamed, and proprietary
closed-source tools to build Cloupe files are disabled.

It is still suffient to generate summary reports and count matrices compatible with
single-cell analysis tools available for 10X Genomics and Cell Ranger output format
in Python and R packages.

## Usage

### Generating References

The Cell Ranger modules can be used to generate reference indexes to run UniverSC.
Note that UniverSC requires the Open Source version v3.0.2 of Cell Ranger included
in the nf-core/universc Docker image. The same module parameters can be run provided
that the container is changed in process configurations (modify nextflow.config).

```
process {

...
    withName: CELLRANGER_MKGTF {
        container = "nf-core/universc:1.2.5.1"
    }
    withName: CELLRANGER_MKREF{
       container = "nf-core/universc:1.2.5.1"
    }
...
}
```

This will generate a compatible index for UniverSC using the same version of the
STAR aligner and a permissive software license without and EULA.

### Container settings

The cellranger install directory must have write permissions to run UniverSC.
To run in docker or podman use the `--user root` option in container parameters
and for singularity use the `--writeable` parameter.

These are set as default in universc/main.nf:

```
    container "nf-core/universc:1.2.5.1"
    if (workflow.containerEngine == 'docker'){
        containerOptions = "--privileged"
    }
    if (workflow.containerEngine == 'podman'){
        containerOptions = "--runtime /usr/bin/crun --userns=keep-id --user root --systemd=always"
    }
    if (workflow.containerEngine == 'singularity'){
        containerOptions = "--writable"
    }
```

Select the container engine with `nextflow --profile "docker"` or set the environment variable
as one of the following before running nextflow.

```
export PROFILE="docker"
export PROFILE="podman"
export PROFILE="singularity"
```

Note that due to dependencies installed in a docker image, it is not possible to use conda environments.

## Disclaimer

We are third party developers not affiliated with 10X Genomics or any other vendor of
single-cell technologies. We are releasing this code on an open-source license which calls Cell Ranger
as an external dependency.

## Licensing

This package is provided open-source on a GPL-3 license. This means that you are free to use and
modify this code provided that they also contain this license.

## Updating the package

The tomkellygenetics/universc:<VERSION> container is automatically updated with tomkellygenetics/universc:latest.

A stable release is mirrored at nf-core/universc:1.2.5.1 and will be updated as needed.

To build an updated container use the Dockerfile provided here:

[https://github.com/minoda-lab/universc/blob/master/Dockerfile](https://github.com/minoda-lab/universc/blob/master/Dockerfile)

Note that this uses a custom base image which is built with an open-source implementation of
Cell Ranger v3.0.2 on MIT License and relies of Python 2. The build file can be found here:

[https://github.com/TomKellyGenetics/cellranger_clean/blob/master/Dockerfile](https://github.com/TomKellyGenetics/cellranger_clean/blob/master/Dockerfile)
