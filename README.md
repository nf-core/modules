# ![nf-core/modules](docs/images/nfcore-modules_logo.png)

[![Nextflow](https://img.shields.io/badge/nextflow-DSL2-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

![GitHub Actions Coda Linting](https://github.com/nf-core/modules/workflows/Code%20Linting/badge.svg)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23modules-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/modules)

[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

> THIS REPOSITORY IS UNDER ACTIVE DEVELOPMENT. SYNTAX, ORGANISATION AND LAYOUT MAY CHANGE WITHOUT NOTICE!

A repository for hosting [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) module files containing tool-specific process definitions and their associated documentation.

## Table of contents

- [Using existing modules](#using-existing-modules)
- [Adding a new module file](#adding-a-new-module-file)
    - [Module template](#module-template)
    - [Guidelines](#guidelines)
    - [Testing](#testing)
    - [Documentation](#documentation)
    - [Uploading to `nf-core/modules`](#uploading-to-nf-coremodules)
- [Terminology](#terminology)
- [Help](#help)
- [Citation](#citation)

## Using existing modules

The module files hosted in this repository define a set of processes for software tools such as `fastqc`, `bwa`, `samtools` etc. This allows you to share and add common functionality across multiple pipelines in a modular fashion.

We have written a helper command in the `nf-core/tools` package that uses the GitHub API to obtain the relevant information for the module files present in the [`software/`](software/) directory of this repository. This includes using `git` commit hashes to track changes for reproducibility purposes, and to download and install all of the relevant module files.

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

4. Import the module in your Nextflow script:

    ```nextflow
    #!/usr/bin/env nextflow

    nextflow.enable.dsl = 2

    include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: [:] )
    ```

5. We have plans to add other utility commands to help developers install and maintain modules downloaded from this repository so watch this space!

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

## Adding a new module file

> **NB:** The definition and standards for module files are still under discussion
but we are now gladly accepting submissions :)

If you decide to upload a module to `nf-core/modules` then this will
ensure that it will become available to all nf-core pipelines,
and to everyone within the Nextflow community! See
[`software/`](software)
for examples.

### Module template

We have added a directory called [`software/SOFTWARE/TOOL/`](software/SOFTWARE/TOOL/) that serves as a template with which to create your own module submission. Where applicable, we have added extensive `TODO` statements to the files in this directory for general information, to help guide you as to where to make the appropriate changes, and how to make them. If in doubt, have a look at how we have done things for other modules.

```console
.
├── software
│   ├── SOFTWARE
│   │   └── TOOL
│   │       ├── functions.nf          ## Utility functions imported in main module script
│   │       ├── main.nf               ## Main module script
│   │       ├── meta.yml              ## Documentation for module, input, output, params, author
│   │       └── test
│   │           ├── input             ## Soft-link input test data from "tests/"
│   │           ├── main.nf           ## Minimal workflow to test module
│   │           ├── nextflow.config   ## Minimal config to test module
│   │           └── output            ## Upload output files from test for unit testing
```

### Guidelines

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

#### General

- Software that can be piped together SHOULD be added to separate module files
unless there is a run-time, storage advantage in implementing in this way. For example,
using a combination of `bwa` and `samtools` to output a BAM file instead of a SAM file:

    ```bash
    bwa mem | samtools view -B -T ref.fasta
    ```

- Where applicable, the usage and generation of compressed files SHOULD be enforced as input and output, respectively:
    - `*.fastq.gz` and NOT `*.fastq`
    - `*.bam` and NOT `*.sam`

- Where applicable, each module command MUST emit a file `<SOFTWARE>.version.txt` containing a single line with the software's version in the format `<VERSION_NUMBER>` or `0.7.17` e.g.

    ```bash
    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${software}.version.txt
    ```

    If the software is unable to output a version number on the command-line then a variable called `VERSION` can be manually specified to create this file e.g. [homer/annotatepeaks module](https://github.com/nf-core/modules/blob/master/software/homer/annotatepeaks/main.nf).

- The process definition MUST NOT contain a `when` statement.

#### Naming conventions

- The directory structure for the module name must be all lowercase e.g. [`software/bwa/mem/`](software/bwa/mem/). The name of the software (i.e. `bwa`) and tool (i.e. `mem`) MUST be all one word.

- The process name in the module file MUST be all uppercase e.g. `process BWA_MEM {`. The name of the software (i.e. `BWA`) and tool (i.e. `MEM`) MUST be all one word separated by an underscore.

- All parameter names MUST follow the `snake_case` convention.

- All function names MUST follow the `camelCase` convention.

#### Module parameters

- A module file SHOULD only define input and output files as command-line parameters to be executed within the process.

- All other parameters MUST be provided as a string i.e. `options.args` where `options` is a Groovy Map that MUST be provided via the Nextflow `addParams` option when including the module via `include` in the parent workflow.

- If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable e.g. `--threads $task.cpus`.

- Any parameters that need to be evaluated in the context of a particular sample e.g. single-end/paired-end data MUST also be defined within the process.

#### Input/output options

- Named file extensions MUST be emitted for ALL output channels e.g. `path "*.txt", emit: txt`.

- Optional inputs are not currently supported by Nextflow. However, "fake files" MAY be used to work around this issue.

#### Resource requirements

- An appropriate resource `label` MUST be provided for the module as listed in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/%7B%7Bcookiecutter.name_noslash%7D%7D/conf/base.config#L29) e.g. `process_low`, `process_medium` or `process_high`.

- If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable e.g. `--threads $task.cpus`.

#### Software requirements

[BioContainers](https://biocontainers.pro/#/) is a registry of Docker and Singularity containers automatically created from all of the software packages on [Bioconda](https://bioconda.github.io/). Where possible we will use BioContainers to fetch pre-built software containers and Bioconda to install software using Conda.

- Software requirements SHOULD be declared within the module file using the Nextflow `container` directive e.g. go to the [BWA BioContainers webpage](https://biocontainers.pro/#/tools/bwa), click on the `Pacakages and Containers` tab, sort by `Version` and get the portion of the link after the `docker pull` command where `Type` is Docker. You may need to double-check that you are using the latest version of the software because you may find that containers for older versions have been rebuilt more recently.

- If the software is available on Conda it MUST also be defined using the Nextflow `conda` directive. Software MUST be pinned to the channel (i.e. `bioconda`) and version (i.e. `0.7.17`) e.g. `bioconda::bwa=0.7.17`. Pinning the build too is not currently a requirement e.g. `bioconda::bwa=0.7.17=h9402c20_2`.

- If required, multi-tool containers may also be available on BioContainers e.g. [`bwa` and `samtools`](https://biocontainers.pro/#/tools/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40). It is also possible for a multi-tool container to be built and added to BioContainers by submitting a pull request on their [`multi-package-containers`](https://github.com/BioContainers/multi-package-containers) repository.

- If the software is not available on Bioconda a `Dockerfile` MUST be provided within the module directory. We will use GitHub Actions to auto-build the containers on the [GitHub Packages registry](https://github.com/features/packages).

#### Publishing results

The [Nextflow `publishDir`](https://www.nextflow.io/docs/latest/process.html#publishdir) definition is currently quite limited in terms of parameter/option evaluation. To overcome this, the publishing logic we have implemented for use with DSL2 modules attempts to minimise changing the `publishDir` directive (default: `params.outdir`) in favour of constructing and appending the appropriate output directory paths via the `saveAs:` statement e.g.

```nextflow
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
```

The `saveFiles` function can be found in the [`functions.nf`](software/fastqc/functions.nf) file of utility functions that will be copied into all module directories. It uses the various publishing `options` specified as input to the module to construct and append the relevant output path to `params.outdir`.

We also use a standardised parameter called `params.publish_dir_mode` that can be used to alter the file publishing method (default: `copy`).

### Testing

- All test data for `nf-core/modules` MUST be added to [`tests/data/`](tests/data/) and organised by filename extension.

- Test files MUST be kept as tiny as possible.

- Every module MUST be tested by adding a test workflow with a toy dataset in the [`test/`](software/fastqc/test) directory of the module.

- Generic files from [`tests/data/`](tests/data/) MUST be reused by symlinking them into the [`test/input/`](software/fastqc/test/input/) directory of the module.

- Any outputs produced by the test workflow MUST be placed in a folder called [`test/output/`](software/fastqc/test/output/) so that they can be used for unit testing.

- If the appropriate test data doesn't exist for your module then it MUST be added to [`tests/data/`](tests/data/).

- A GitHub Actions workflow file MUST be added to [`.github/workflows/`](.github/workflows/) e.g. [`.github/workflows/fastqc.yml`](.github/workflows/fastqc.yml).

### Documentation

- A module MUST be documented in the [`meta.yml`](software/fastqc/meta.yml) file. It MUST document `params`, `input` and `output`. `input` and `output` MUST be a nested list.

### Uploading to `nf-core/modules`

[Fork](https://help.github.com/articles/fork-a-repo/) the `nf-core/modules` repository to your own GitHub account. Within the local clone of your fork add the module file to the [`software/`](software) directory.

Commit and push these changes to your local clone on GitHub, and then [create a pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) on the `nf-core/modules` GitHub repo with the appropriate information.

We will be notified automatically when you have created your pull request, and providing that everything adheres to nf-core guidelines we will endeavour to approve your pull request as soon as possible.

## Terminology

The features offered by Nextflow DSL2 can be used in various ways depending on the granularity with which you would like to write pipelines. Please see the listing below for the hierarchy and associated terminology we have decided to use when referring to DSL2 components:

- *Module*: A `process` that can be used within different pipelines and is as atomic as possible i.e. cannot be split into another module. An example of this would be a module file containing the process definition for a single tool such as `FastQC`. At present, this repository has been created to only host atomic module files that should be added to the [`software/`](software/) directory along with the required documentation and tests.

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


<!---

- The module MUST accept a parameter `params.publish_results` accepting at least
    - `"none"`, to publish no files at all,
    - a glob pattern which is initalized to a sensible default value.

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

-->
