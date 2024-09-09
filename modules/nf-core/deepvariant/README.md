# Conda is not supported at the moment

The [bioconda](https://bioconda.github.io/recipes/deepvariant/README.html) recipe is not fully working as expected.

See https://github.com/bioconda/bioconda-recipes/issues/30310 and https://github.com/nf-core/modules/issues/1754 for more information.

Hence, we are using the docker container provided by the authors of the tool:

- [google/deepvariant](https://hub.docker.com/r/google/deepvariant)

This image is mirrored on the [nf-core quay.io](https://quay.io/repository/nf-core/deepvariant) for convenience.
