# DeepVariant subworkflow

Usage: the input channel should contain tuples of three elements: `meta`, an alignment file in bam or
cram format, and a corresponding index.

It is very important that `meta` is unique for all the input elements, because the pipeline does a join on `meta`.

Please note the important configuration items listed in the `deepvariant` module's README file. It is required to specify the input channels for `DEEPVARIANT_MAKEEXAMPLES` and the model to run for `DEEPVARIANT_CALLVARIANTS`.

