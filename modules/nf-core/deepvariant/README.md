# DeepVariant module / subworkflow options

The DeepVariant tool can be run using the `deepvariant/rundeepvariant` subcommand, or the subworkflow `deepvariant`, which calls the subcommands `makeexamples`, `callvariants` and `postprocessvariants`.  The subcommand `rundeepvariant` is simpler, but the subworkflow may be useful if you want to run `callvariants` on GPU.


# Conda is not supported at the moment

The [bioconda](https://bioconda.github.io/recipes/deepvariant/README.html) recipe is not fully working as expected.

See https://github.com/bioconda/bioconda-recipes/issues/30310 and https://github.com/nf-core/modules/issues/1754 for more information.

Hence, we are using the docker container provided by the authors of the tool:

- [google/deepvariant](https://hub.docker.com/r/google/deepvariant)

This image is mirrored on the [nf-core quay.io](https://quay.io/repository/nf-core/deepvariant) for convenience.


# DeepVariant subworkflow

You can use the subworkflow `nf-core/deepvariant`, which integrates the three
processes to perform variant calling with common file formats.

These module subcommands incorporate the individual steps of the DeepVariant pipeline:

    * makeexamples: Converts the input alignment file to a tfrecord format suitable for the deep learning model
    * callvariants: Call variants based on input tfrecords. The output is also in
    tfrecord format, and needs postprocessing to convert it to vcf.
    * postprocessvariants: Convert variant calls from callvariants to VCF, and
    also create GVCF files based on genomic information from makeexamples.

