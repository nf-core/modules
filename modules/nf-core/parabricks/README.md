# NVIDIA Parabricks Modules

These nf-core modules implement functionality of the NVIDIA Parabricks software for GPU-accelerated genomics tasks. Parabricks covers the functionality of alignment (Ex. `bwa mem`, `bwa meth` and `STAR`), variant and fusion calling (Ex. `deepvariant`, `mutectcaller`, and `starfusion`), and other common genomics tasks. Please see the [documentation](https://docs.nvidia.com/clara/parabricks/latest/index.html) for additional details.

## General considerations

Parabricks is available only through a docker container: `nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1`. The main entrypoint to the Parabricks is the command line program `pbrun`, which calls several sub-programs within the container. The Parabricks authors combined several common patterns into a single tool, such as `fq2bam` performing alignment, sorting, mark duplicates, and base quality score recalibration, all within a single command. Generally, the tools can be used for only a subset of the entire functionality as well.

## Hardware Requirements

Parabricks tools must be run with at least one NVIDIA GPU with at least 16GB vRAM, and a usually require a large amount of resources (at least 8 threads and 30GB RAM). Please see the [installation requirements](https://docs.nvidia.com/clara/parabricks/latest/gettingstarted/installationrequirements.html) for more information. It is recommended to use fast local storage whenever possible to ensure the system is not bottlenecked by I/O.

To give docker or singularity access to GPUs present on the host system, add this line to the configuration file: `docker.runOptions = "--gpus all"` or `singularity.runOptions = "--nv"`.

## License

As of version 4.0, Parabricks is available for free for general use. The license from NVIDIA states:

```
A license is no longer required to use Clara Parabricks. The container works out of the box once downloaded.

Clara Parabricks is free for

    * Academic research, or
    * Research by non-profit institutions, or
    * For development, test and evaluation purposes without use in production.

Users who would like to have Enterprise Support for Clara Parabricks can purchase NVIDIA AI Enterprise licenses, which provides full-stack support. To learn more about NVIDIA AI Enterprise, please visit https://www.nvidia.com/en-us/data-center/products/ai-enterprise/.
```

## Notes on Specific tools

### fq2bam

The `fq2bam` module performs alignment, sorting, (optional) marking of duplicates, and (optional) base quality score recalibration (BQSR). There is no option to control the number of threads used with this tool - all available threads on the system are used by default.

Alignment and coordinate sorting are always performed. Duplicate marking can be performed by passing the option `markdups=true`. Duplicate marking and BQSR can be performed by passing the options `markdups=true` and `known_sites=$KNOWN_SITES_FILE`.

Please see the `fq2bam/meta.yml` file for a detailed list of required and optional inputs and outputs.

For additional considerations, including information about how readgroups are added to the resulting bam files, see the [tool documentation](https://docs.nvidia.com/clara/parabricks/latest/Documentation/ToolDocs/man_fq2bam.html).

## rnafq2bam

The `rnafq2bam` module is based on STAR version `2.7.2a`. Therefore the genome lib directory required as input for this module must also be generated using this version of STAR. For convenience, a module with this version of STAR is included in this directory under `parabricks/stargenomegenerate`.

## starfusion

The `starfusion` module is based on starfusion 1.7.0. Therefore the genome lib directory required as input for this module must also be generated using this version of starfusion. For convenience, a module with this version of starfusion is included in this directory under `parabricks/starfusion_build`.

## compatible_with.yaml

The `compatible_with.yaml` is provided as an optional output to the stub section to make the compatible CPU version available to the end user. This section is not given for the subtools `applybqsr`, `fq2bammeth`, `genotypegvcf`, `rnafq2bam`, or `starfusion`.

For the full list of compatible versions, check the [Parabricks documentation](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/outputaccuracyandcompatiblecpusoftwareversions.html#).

## Notes on Testing

The following Parabricks modules require testing on a `g4dn.2xlarge` instead of the default `g4dn.xlarge` due to higher system memory requirements. 

* rnafq2bam
* starfusion

In [nf-test-gpu.yml](../../../.github/workflows/nf-test-gpu.yml), update the following line: 

```
nf-test-gpu:
    runs-on: "runs-on=${{ github.run_id }}/family=g4dn.2xlarge/image=ubuntu24-gpu-x64"
```

Change it back before merging into master:

```
nf-test-gpu:
    runs-on: "runs-on=${{ github.run_id }}/family=g4dn.xlarge/image=ubuntu24-gpu-x64"
```