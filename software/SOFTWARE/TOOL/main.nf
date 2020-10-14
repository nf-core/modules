// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// TODO nf-core: All of these TODO statements can be deleted after the relevant changes have been made.
// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/software
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join

// TODO nf-core: The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in RFC 2119 (https://tools.ietf.org/html/rfc2119).
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided as a string i.e. "options.args"
//               where "params.options" is a Groovy Map that MUST be provided via the addParams section of the including workflow.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. bwa mem | samtools view -B -T ref.fasta to output BAM instead of SAM.
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, "fake files" MAY be used to work around this issue.

params.options = [:]
def options    = initOptions(params.options)

// TODO nf-core: Process name MUST be all uppercase,
//               "SOFTWARE" and (ideally) "TOOL" MUST be all one word separated by an "_".
process SOFTWARE_TOOL {
    // TODO nf-core: If a meta map of sample information is NOT provided in "input:" section
    //               change tag value to another appropriate input value e.g. tag "$fasta"
    tag "$meta.id"
    // TODO nf-core: Provide appropriate resource label for process as listed in the nf-core pipeline template below:
    //               https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/%7B%7Bcookiecutter.name_noslash%7D%7D/conf/base.config#L29
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        // TODO nf-core: If a meta map of sample information is NOT provided in "input:" section
        //               change "publish_id:meta.id" to initialise an empty string e.g. "publish_id:''".
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    // TODO nf-core: List required Conda packages.
    //               Software MUST be pinned to channel (i.e. "bioconda") and version (i.e. "1.10") as in the example below.
    //               Pinning the build too e.g. "bioconda::samtools=1.10=h9402c20_2" is not currently a requirement.
    conda (params.enable_conda ? "bioconda::samtools=1.10" : null)

    // TODO nf-core: Fetch "docker pull" address for latest BioContainer image of software: e.g. https://biocontainers.pro/#/tools/samtools
    //               Click on the Pacakages and Containers tab, sort by Version and get the portion of the link after the docker pull command where Type is Docker.
    //               You may need to double-check that you are using the latest version of the software because you may find that containers for older versions have been rebuilt more recently.
    //               If required, multi-tool containers may also be available and are usually named to start with "mulled".
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"
    } else {
        container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    }

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/software/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(reads)
    
    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    // TODO nf-core: If meta is provided in "input:" section then it MUST be added to ALL output channels (except version)
    tuple val(meta), path("*.bam"), emit: bam
    // TODO nf-core: List additional required output channels/values here
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    // TODO nf-core: If a meta map of sample information is NOT provided in "input:" section delete the line below
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/software/homer/annotatepeaks/main.nf
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$options.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    software tool \\
        $options.args \\
        --threads $task.cpus \\
        $reads \\
        > ${prefix}.bam

    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}
