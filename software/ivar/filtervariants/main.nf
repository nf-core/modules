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

//Usage: ivar filtervariants -p <prefix> replicate-one.tsv replicate-two.tsv ... OR ivar filtervariants -p <prefix> -f <text file with one variant file per line> 
// Input: Variant tsv files for each replicate/sample
// Input Options    Description
//            -t    Minimum fration of files required to contain the same variant. Specify value within [0,1]. (Default: 1)
//            -f    A text file with one variant file per line.
// Output Options   Description
//            -p    (Required) Prefix for the output filtered tsv file



params.options = [:]
def options    = initOptions(params.options)

process IVAR_FILTERVARIANTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::ivar=1.3" : null)

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ivar:1.3--h089eab3_1"
    } else {
        container "quay.io/biocontainers/ivar:1.3--h089eab3_1"
    }

    input:
    tuple val(meta), path(variant_tsv)

    output:
    tuple val(meta), path("*.filtered.tsv"), emit: filtered_tsv
    path "*.version.txt"                   , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    // is the input a file of filenames (one per line)?
    def fofn     = options.fofn ? "-f " : ""
    """
    ivar filtervariants \\
        $options.args \\
        -p $prefix \\
        ${fofn}${variant_tsv}

    ivar version | head -n1 2>&1 | sed 's/^.*iVar version //g' > ${software}.version.txt
    """
}
