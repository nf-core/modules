# ![nf-core/modules](docs/images/nfcore-modules_logo.png)

> DSL2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT. SYNTAX, ORGANISATION AND LAYOUT OF THIS REPOSITORY MAY CHANGE IN THE NEAR FUTURE!

A repository for hosting nextflow [`DSL2`](https://www.nextflow.io/docs/edge/dsl2.htmlhttps://www.nextflow.io/docs/edge/dsl2.html) module files and subworkflows containing tool-specific process definitions and their associated documentation.

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

The DSLv2 features for Nextflow are new for everyone and not just beginners. We discussed some terminology terms to discuss things related to modules & subworkflows more appropriately:

* *Module*: A `process`that can be used between different pipelines, which is atomic (i.e. can/should not be divided further).
* *Subworkflow*: A combined set of modules, that combine a logical step in a pipeline using multiple modules together. Example: A preprocessing subworkflow for FastQ input, that could be (re-) used across multiple pipelines to QC input FastQ files.
* *Workflow*: What DSLv1 users would consider a pipeline, e.g. from input to potentially complex output. Can either consist of individual modules, a large monolithic script as in DSLv1 or a combination of subworkflows (or any combination of these three). 

## Organization of repository

* Modules should end up in the subdirectory `tools`
* Subworkflows should be kept with individual pipelines, and not end up here in the modules repository. 

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
