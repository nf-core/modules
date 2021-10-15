// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/software
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join

// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided as a string i.e. "options.args"
//               where "params.options" is a Groovy Map that MUST be provided via the addParams section of the including workflow.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

params.options = [:]
options        = initOptions(params.options)

process IMPUTEME_VCFTOPRS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "YOUR-TOOL-HERE" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "quay.io/lassefolkersen/imputeme:v1.0.6"
    } else {
        container "quay.io/lassefolkersen/imputeme:v1.0.6"
    }

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/software/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(vcf)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.json"), emit: json
    // TODO nf-core: List additional required output channels/values here
    path "*.version.txt"          , emit: version


    containerOptions "\
        --mount 'type=tmpfs,source=,target=/home/ubuntu/logs' \
        --mount 'type=tmpfs,source=,target=/home/ubuntu/misc_files' \
        --mount 'type=volume,source=,target=/home/ubuntu/configuration' \
        --mount 'type=tmpfs,source=,target=/home/ubuntu/data' \
        --mount 'type=volume,source=,target=/home/ubuntu/programs' \
        --mount 'type=volume,source=,target=/home/ubuntu/prs_dir' \
        --mount 'type=tmpfs,source=,target=/home/ubuntu/imputations' \
        --mount 'type=tmpfs,source=,target=/home/ubuntu/vcfs' \
        --mount 'type=volume,source=,target=/home/ubuntu/srv'"


    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/software/homer/annotatepeaks/main.nf
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$options.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    #!/usr/bin/env Rscript
    dir.create("~/logs/submission")
    source("/home/ubuntu/srv/impute-me/functions.R")
    prepare_individual_genome('$vcf')
    library(jsonlite)
    list(test0=list.files("/home/ubuntu/test/ubuntu"),home=list.files("/home/ubuntu"),root=list.files("/"),test=list.files("test"),test2=list.files("/home/ubuntu/test"),test3=list.files("~/srv"),test4=list.files("~/configuration"))->d
    filename<-"test.json"
    JSON<-toJSON(d,digits=NA)
    f1<-file(filename,"w")
    writeLines(JSON,f1)
    close(f1)


    version_file_path="${software}.version.txt"
    f2<-file(version_file_path,"w")
    writeLines("v1.0.6",f2)
    close(f2)
    #echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}
