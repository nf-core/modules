# ![nf-core/modules](docs/images/nfcore-modules_logo.png)

![GitHub Actions Coda Linting](https://github.com/nf-core/modules/workflows/Code%20Linting/badge.svg)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23modules-4A154B?logo=slack)](https://nfcore.slack.com/channels/modules)

> THIS REPOSITORY IS UNDER ACTIVE DEVELOPMENT. SYNTAX, ORGANISATION AND LAYOUT MAY CHANGE WITHOUT NOTICE!

A repository for hosting Nextflow [`DSL2`](https://www.nextflow.io/docs/latest/dsl2.html) module files (see [Terminology](#terminology)) containing tool-specific process definitions and their associated documentation.

## Table of contents

- [Using existing modules](#using-existing-modules)
    - [Configuration and parameters](#configuration-and-parameters)
    - [Offline usage](#offline-usage)
- [Adding a new module file](#adding-a-new-module-file)
    - [Testing](#testing)
    - [Documentation](#documentation)
    - [Uploading to `nf-core/modules`](#uploading-to-nf-coremodules)
- [Terminology](#terminology)
- [Help](#help)
- [Citation](#citation)

## Using existing modules

We have written a helper command in the `nf-core/tools` package that allows you to install any module present in the `software/` directory of this repository:

1. [Install](https://github.com/nf-core/tools#installation) the latest version of `nf-core/tools` (`>=1.10.2`)
2. List the available modules:

    ```console
    $ nf-core modules list

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 1.10.2



    INFO      Modules available from nf-core/modules (master):                                                                                                                  modules.py:51

    bwa/index
    bwa/mem
    deeptools/computematrix
    deeptools/plotfingerprint
    deeptools/plotheatmap
    deeptools/plotprofile
    fastqc
    ..truncated..
    ```

3. Install the module in your pipeline directory:

    ```console
    $ nf-core modules install . fastqc

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 1.10.2



    INFO      Installing fastqc                                                                                                                                                 modules.py:62
    INFO      Downloaded 3 files to ./modules/nf-core/software/fastqc                                                                                                           modules.py:97
    ```

4. We have plans to add other utility commands to help developers install and maintain modules downloaded from this repository so watch this space!

    ```console
    $ nf-core modules --help

    ...truncated...

    Commands:
      list     List available software modules.
      install  Add a DSL2 software wrapper module to a pipeline.
      update   Update one or all software wrapper modules.             (NOT YET IMPLEMENTED)
      remove   Remove a software wrapper from a pipeline.              (NOT YET IMPLEMENTED)
      check    Check that imported module code has not been modified.  (NOT YET IMPLEMENTED)
```

### Configuration and parameters

The module files hosted in this repository define a set of processes for software tools such as `fastqc`, `trimgalore`, `bwa` etc. This allows you to share and add common functionality across multiple pipelines in a modular fashion.

> The definition and standards for module files are still under discussion amongst the community but hopefully, a description should be added here soon!

### Offline usage

If you want to use an existing module file available in `nf-core/modules`, and you're running on a system that has no internet connection, you'll need to download the repository (e.g. `git clone https://github.com/nf-core/modules.git`) and place it in a location that is visible to the file system on which you are running the pipeline. Then run the pipeline by creating a custom config file called e.g. `custom_module.conf` containing the following information:

```bash
include /path/to/downloaded/modules/directory/
```

Then you can run the pipeline by directly passing the additional config file with the `-c` parameter:

```bash
nextflow run /path/to/pipeline/ -c /path/to/custom_module.conf
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs + modules in one go for you, to make this process easier.

## Adding a new module file

If you decide to upload your module file to `nf-core/modules` then this will
ensure that it will be automatically downloaded, and available at run-time to
all nf-core pipelines, and to everyone within the Nextflow community! See
[`nf-core/modules/software`](https://github.com/nf-core/modules/tree/master/software)
for examples.

**The definition and standards for module files are still under discussion
amongst the community. Currently the following points have been agreed on:**

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

### Defining inputs, outputs and parameters
- A module file SHOULD only define inputs and outputs as parameters. Additionally,
    - it MUST define threads or resources where required for a particular process using `task.cpus`
    - ~~it MUST be possible to pass additional parameters to the tool as a command line string via the `params.<MODULE>_args` parameter.~~
    - it MUST be possible to pass additional parameters as a [nextflow Map](https://www.nextflow.io/docs/latest/script.html#maps) through an additional input channel `val(options)` [Details require discussion].
    - All NGS modules MUST accept a triplet [name, single_end, reads] as input. The single-end boolean values MUST be specified through the input channel and not inferred from the data e.g. [here](https://github.com/nf-core/tools/blob/028a9b3f9d1ad044e879a1de13d3c3a25a06b9a7/nf_core/pipeline-template/%7B%7Bcookiecutter.name_noslash%7D%7D/modules/nf-core/fastqc.nf#L13).
- Process names MUST be all uppercase.
- Each process MUST emit a file `<TOOL>.version.txt` containing a single line with the software's version in the format `v<VERSION_NUMBER>`.
- All outputs MUST be named using `emit`.
- A Process MUST NOT contain a `when` statement.
- Optional inputs need development on the nextflow side. In the meanwhile, "fake files" MAY be used to work around this issue.

### Atomicity
- Software that can be piped together SHOULD be added to separate module files unless there is an run-time, storage advantage in implementing in this way e.g. `bwa mem | samtools view -C -T ref.fasta` to output CRAM instead of SAM.

### Resource requirements
- Each module MUST define a label `process_low`, `process_medium` or `process_high` to declare resource requirements. (*These flags will be ignored outside of nf-core and the pipeline developer is free to define adequate resource requirements*)

### Publishing results
- The module MUST accept the parameters `params.out_dir` and `params.publish_dir` and MUST publish results into `${params.out_dir}/${params.publish_dir}`.
- The `publishDirMode` MUST be configurable via `params.publish_dir_mode`
- The module MUST accept a parameter `params.publish_results` accepting at least
    - `"none"`, to publish no files at all,
    - a glob pattern which is initalized to a sensible default value.

    It MAY accept `"logs"` to publish relevant log files, or other flags, if applicable.

- To ensure consistent naming, files SHOULD be renamed according to the `$name` variable before returning them.

### Testing
- Every module MUST be tested by adding a test workflow with a toy dataset.
- Test data MUST be stored within this repo. It is RECOMMENDED to re-use generic files from `tests/data` by symlinking them into the test directory of the module. Specific files MUST be added to the test-directory directly. Test files MUST be kept as tiny as possible.

### Software requirements
- Software requirements SHOULD be declared in a conda `environment.yml` file, including exact version numbers. Additionally, there MUST be a `Dockerfile` that containerizes the environment, or packages the software if conda is not available.
- Docker containers MUST BE identified by their `sha256(Dockerfile + environment.yml)`.
- Each module must have it's own `Dockerfile` and `environment.yml` file
    - Care should be taken to maintain identical files for subcommands that use the same software. Then the hash tag will be the same and they will be implicitly re-used across subcommands.

### File formats
- Wherever possible, [CRAM](https://en.wikipedia.org/wiki/CRAM_(file_format)) files SHOULD be used over BAM files.
- Wherever possible, FASTQ files SHOULD be compressed using gzip.

### Documentation
- A module MUST be documented in the `meta.yml` file. It MUST document `params`, `input` and `output`. `input` and `output` MUST be a nested list. [Exact detail need to be elaborated. ]

### Uploading to `nf-core/modules`

[Fork](https://help.github.com/articles/fork-a-repo/) the `nf-core/modules` repository to your own GitHub account. Within the local clone of your fork add the module file to the [`nf-core/modules/software`](https://github.com/nf-core/modules/tree/master/software) directory. Please keep the naming consistent between the module and documentation files e.g. `bwa.nf` and `bwa.md`, respectively.

Commit and push these changes to your local clone on GitHub, and then [create a pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) on `nf-core/modules` GitHub repo with the appropriate information.

We will be notified automatically when you have created your pull request, and providing that everything adheres to nf-core guidelines we will endeavour to approve your pull request as soon as possible.

## Terminology

The features offered by Nextflow DSL2 can be used in various ways depending on the granularity with which you would like to write pipelines. Please see the listing below for the hierarchy and associated terminology we have decided to use when referring to DSL2 components:

- *Module*: A `process` that can be used within different pipelines and is as atomic as possible i.e. cannot be split into another module. An example of this would be a module file containing the process definition for a single tool such as `FastQC`. At present, this repository has been created to only host atomic module files that should be added to the `software/` directory along with the required documentation and tests.
- *Sub-workflow*: A chain of multiple modules that offer a higher-level of functionality within the context of a pipeline. For example, a sub-workflow to run multiple QC tools with FastQ files as input. Sub-workflows should be shipped with the pipeline implementation and if required they should be shared amongst different pipelines directly from there. As it stands, this repository will not host sub-workflows although this may change in the future since well-written sub-workflows will be the most powerful aspect of DSL2.
- *Workflow*: What DSL1 users would consider an end-to-end pipeline. For example, from one or more inputs to a series of outputs. This can either be implemented using a large monolithic script as with DSL1, or by using a combination of DSL2 individual modules and sub-workflows.

## Help

For further information or help, don't hesitate to get in touch on [Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

If you use the module files in this repository for your analysis please you can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
