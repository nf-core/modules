# DeepVariant subworkflow

Usage: the input channel should contain tuples of three elements: `meta`, an alignment file in bam or
cram format, and a corresponding index.

It is very important that the input channel's `meta` is unique for all the input elements, because the subworkflow does a join on `meta`.

Please note the important configuration items listed in the `deepvariant` module's README file. It is required to use the configuration to specify the input "channels" (data types to extract from bam file) for `DEEPVARIANT_MAKEEXAMPLES`, and the model to run for `DEEPVARIANT_CALLVARIANTS`. The correct arguments for a specific model (data type) can be determined by manually using the `run_deepvariant` command from the Docker / Singularity image with the `--dry_run` option.
