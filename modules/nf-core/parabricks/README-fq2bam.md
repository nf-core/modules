# Parabricks fq2bam

`fq2bam` performs alignment, sorting, (optional) marking of duplicates, and (optional) base quality score recalibration (BQSR). There is no option to control the number of threads used with this tool - all avilable threads on the system are used by default.

Alignment and coordinate sorting are always performed. Duplicate marking can be performed by passing the option `markdups=true`. Duplicate marking and BQSR can be performed by passing the options `markdups=true` and `known_sites=$KNOWN_SITES_FILE`.

Please see the `fq2bam/meta.yml` file for a detailed list of required and optional inputs and outputs.

For additional considerations, including information about how readgroups are added to the resulting bam files, see the [tool documentation](https://docs.nvidia.com/clara/parabricks/4.0.1/Documentation/ToolDocs/man_fq2bam.html).
