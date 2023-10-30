# NVIDIA Clara Parabricks modules

These nf-core modules implement functionality of the NVIDIA Clara Parabricks programs for GPU-accelerated genomics tasks. Parabricks covers the functionality of alignment (replicating `bwa mem`), variant calling (replicating `gatk4` and `deepvariant`), and other common genomics tasks. Please see the [documentation](https://docs.nvidia.com/clara/parabricks/4.0.1/index.html) for additional details.

## General considerations

Parabricks is available only through a docker container: `nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1`. The main entrypoint to the Parabricks is the command line program `pbrun`, which calls several sub-programs within the container. The Prabricks authors combined several common patterns into a single tool, such as `fq2bam` performing alignment, sorting, mark duplicates, and base quality score recalibration, all within a single command. Generally, the tools can be used for only a subset of the entire functionality as well.

Parabricks tools must be run with at least one NVIDIA GPU with >16GB vRAM, and a usually require a large amount of resources (at least 8 threads and 30GB RAM). Please see the [resource requirements](https://docs.nvidia.com/clara/parabricks/4.0.1/GettingStarted.html) for more information. It is recommended to use fast local storage whenever possible to ensure the system is not bottlenecked by I/O.

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

Please see the README file for each individual tool for more information.
