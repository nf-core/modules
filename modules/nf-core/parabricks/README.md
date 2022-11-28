# NVIDIA Clara Parabricks modules

These nf-core modules implement functionality of the NVIDIA Clara Parabricks programs for GPU-accelerated genomics tasks. Parabricks covers the functionality of alignment (replicating `bwa mem`), variant calling (replicating `gatk4` and `deepvariant`), and other common genomics tasks. Please see the [documentation](https://docs.nvidia.com/clara/parabricks/4.0.0/index.html) for additional details.

## General considerations

Parabricks is available only through a docker container: `nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1`. The main entrypoint to the Parabricks is the command line program `pbrun`, which calls several sub-programs within the container. The Prabricks authors combined several common patterns into a single tool, such as `fq2bam` performing alignment, sorting, mark duplicates, and base quality score recalibration, all within a single command. Generally, the tools can be used for only a subset of the entire functionality as well.

Parabricks tools must be run with at least one NVIDIA GPU with >16GB vRAM, and a usually require a large amount of resources (at least 8 threads and 30GB RAM). Please see the [resource requirements](https://docs.nvidia.com/clara/parabricks/4.0.0/GettingStarted.html) for more information. It is recommended to use fast local storage whenever possible to ensure the system is not bottlenecked by I/O.

To give docker or singularity access to GPUs present on the host system, add this line to the configuration file: `docker.runOptions = "--gpus all"` or `singularity.runOptions = "--nv"`.

## License

As of version 4.0, Parabricks is available for general use. The license from NVIDIA states:

```
A license is no longer required to use Clara Parabricks. The container works out of the box once downloaded.

Clara Parabricks is free for

    * Academic research, or
    * Research by non-profit institutions, or
    * For development, test and evaluation purposes without use in production.

Users who would like to have Enterprise Support for Clara Parabricks can purchase NVIDIA AI Enterprise licenses, which provides full-stack support. To learn more about NVIDIA AI Enterprise, please visit https://www.nvidia.com/en-us/data-center/products/ai-enterprise/.
```

## Specific tools

### fq2bam

`fq2bam` performs alignment, sorting, (optional) marking of duplicates, and (optional) base quality score recalibration (BQSR). There is no option to control the numebr of threads used with this tool - all avilable threads on the system are used by default.

Alignment and coordinate sorting are always performed. Duplicate marking can be performed by passing the option `markdups=true`. Duplicate marking and BQSR can be performed by passing the options `markdups=true` and `known_sites=$KNOWN_SITES_FILE`.

Please see the `fq2bam/meta.yml` file for a detailed list of required and optional inputs and outputs.

For additional considerations, including information about how readgroups are added to the resulting bam files, see the [tool documentation](https://docs.nvidia.com/clara/parabricks/4.0.0/Documentation/ToolDocs/man_fq2bam.html).

### applybqsr

TBD.

### mutectcaller

TBD.

### haplotypecaller

TBD.

### deepvariant

TBD.
