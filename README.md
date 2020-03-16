# ![nf-core/modules](docs/images/nfcore-modules_logo.png)

> DSL2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT. SYNTAX, ORGANISATION AND LAYOUT OF THIS REPOSITORY MAY CHANGE IN THE NEAR FUTURE!

A repository for hosting nextflow [`DSL2`](https://www.nextflow.io/docs/edge/dsl2.htmlhttps://www.nextflow.io/docs/edge/dsl2.html) module files containing tool-specific process definitions and their associated documentation.

## Table of contents
* [Using existing modules](#using-existing-modules)
    * [Configuration and parameters](#configuration-and-parameters)
    * [Offline usage](#offline-usage)
* [Adding a new module file](#adding-a-new-module-file)
    * [Testing](#testing)
    * [Documentation](#documentation)
    * [Uploading to `nf-core/modules`](#uploading-to-nf-coremodules)
* [Help](#help)

## Terminology

The features offered by Nextflow DSL 2 can be used in various ways depending on the granularity with which you would like to write pipelines. Please see the listing below for the hierarchy and associated terminology we have decided to use when referring to DSL 2 components:

* *Module*: A `process`that can be used within different pipelines and is as atomic as possible i.e. cannot be split into another module. An example of this would be a module file containing the process definition for a single tool such as `FastQC`.
* *Sub-workflow*: A chain of multiple modules that offer a higher-level of functionality within the context of a pipeline.  For example, a sub-workflow to run multiple QC tools with FastQ files as input.
* *Workflow*: What DSL 1 users would consider an end-to-end pipeline. For example, from one or more inputs to a series of outputs. This can either be implemented using a large monolithic script as with DSL 1, or by using a combination of DSL 2 individual modules and sub-workflows. 

## Organization of repository

* This repository has been created to only host atomic module files that should be added to the `tools` sub-directory along with the required documentation, software and tests.
* Sub-workflows should be shipped with the pipeline implementation and if required they should be shared amongst different pipelines directly from there. As it stands, this repository will not host sub-workflows.

## Using existing modules

The Nextflow [`include`](https://www.nextflow.io/docs/edge/dsl2.html#modules-include) statement can be used within your pipelines in order to load module files that you have available locally.

You should be able to get a good idea as to how other people are using module files by looking at pipelines available in nf-core e.g. [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq/pull/162)

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

If you decide to upload your module file to `nf-core/modules` then this will ensure that it will be automatically downloaded, and available at run-time to all nf-core pipelines, and to everyone within the Nextflow community! See [`nf-core/modules/nf`](https://github.com/nf-core/modules/tree/master/nf) for examples.

> The definition and standards for module files are still under discussion amongst the community but hopefully, a description should be added here soon!

### Testing

If you want to add a new module config file to `nf-core/modules` please test that your pipeline of choice runs as expected by using the [`-include`](https://www.nextflow.io/docs/edge/dsl2.html#modules-include) statement with a local version of the module file.

### Documentation

Please add some documentation to the top of the module file in the form of native Nextflow comments. This has to be specified in a particular format as you will be able to see from other examples in the [`nf-core/modules/nf`](https://github.com/nf-core/modules/tree/master/nf) directory.

### Uploading to `nf-core/modules`

[Fork](https://help.github.com/articles/fork-a-repo/) the `nf-core/modules` repository to your own GitHub account. Within the local clone of your fork add the module file to the [`nf-core/modules/nf`](https://github.com/nf-core/modules/tree/master/nf) directory. Please keep the naming consistent between the module and documentation files e.g. `bwa.nf` and `bwa.md`, respectively.

Commit and push these changes to your local clone on GitHub, and then [create a pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) on `nf-core/modules` GitHub repo with the appropriate information.

We will be notified automatically when you have created your pull request, and providing that everything adheres to nf-core guidelines we will endeavour to approve your pull request as soon as possible.

## Help

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).
